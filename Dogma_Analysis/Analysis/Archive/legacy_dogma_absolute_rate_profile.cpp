#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

struct Config {
    std::string inputPath;
    std::string outputPrefix;
    double binWidthNs = 1000.0;
    double triggerResetThresholdSeconds = 1.0;
};

struct BinCounts {
    std::uint32_t ncal1 = 0;
    std::uint32_t lStilbene = 0;
    std::uint32_t sStilbene = 0;
};

struct EventState {
    std::array<std::array<std::array<double, 2>, 32>, 2> times{};
    std::array<std::int64_t, 2> triggerNs{};
    std::array<bool, 2> hasTrigger{};
    bool hasAnyData = false;
};

struct AnalysisStats {
    std::uint64_t parsedHeaders = 0;
    std::uint64_t parsedRows = 0;
    std::uint64_t finalizedEvents = 0;
    bool duplicateRestartDetected = false;
    double duplicateRestartPreviousTriggerSeconds = 0.0;
    double duplicateRestartNewTriggerSeconds = 0.0;
    std::uint64_t duplicateRestartHeadersSkipped = 0;
    bool hasRunStartTrigger = false;
    std::int64_t runStartTriggerNs = 0;
    std::uint64_t validNcal1Hits = 0;
    std::uint64_t validLStilbeneHits = 0;
    std::uint64_t validSStilbeneHits = 0;
    std::int64_t minBinIndex = std::numeric_limits<std::int64_t>::max();
    std::int64_t maxBinIndex = std::numeric_limits<std::int64_t>::min();
    std::int64_t minHitAbsoluteNs = std::numeric_limits<std::int64_t>::max();
    std::int64_t maxHitAbsoluteNs = std::numeric_limits<std::int64_t>::min();
};

constexpr int kTargetTdcIndex = 1;
constexpr int kNcal1Channel = 2;
constexpr int kSmallStilbeneChannel = 20;
constexpr int kLargeStilbeneChannel = 22;
constexpr double kMaxAcceptedHitTimeNs = 50000.0;

bool starts_with(std::string_view text, std::string_view prefix) {
    return text.size() >= prefix.size() && text.substr(0, prefix.size()) == prefix;
}

Config parse_args(int argc, char** argv) {
    Config config;
    for (int index = 1; index < argc; ++index) {
        const std::string arg = argv[index];
        if (arg == "--input" && index + 1 < argc) {
            config.inputPath = argv[++index];
        } else if (arg == "--output-prefix" && index + 1 < argc) {
            config.outputPrefix = argv[++index];
        } else if (arg == "--bin-width-ns" && index + 1 < argc) {
            config.binWidthNs = std::stod(argv[++index]);
        } else if (arg == "--trigger-reset-threshold-seconds" && index + 1 < argc) {
            config.triggerResetThresholdSeconds = std::stod(argv[++index]);
        } else if (arg == "--help") {
            std::cout
                << "Usage: dogma_absolute_rate_profile --input <file> --output-prefix <prefix>\n"
                << "  [--bin-width-ns 1000] [--trigger-reset-threshold-seconds 1.0]\n";
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown or incomplete argument: " + arg);
        }
    }

    if (config.inputPath.empty()) {
        throw std::runtime_error("Missing required argument: --input");
    }
    if (config.outputPrefix.empty()) {
        throw std::runtime_error("Missing required argument: --output-prefix");
    }
    if (config.binWidthNs <= 0.0) {
        throw std::runtime_error("bin-width-ns must be positive");
    }
    if (config.triggerResetThresholdSeconds <= 0.0) {
        throw std::runtime_error("trigger-reset-threshold-seconds must be positive");
    }
    return config;
}

std::string format_double(double value) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(6) << value;
    return stream.str();
}

bool parse_header(const std::string& line, int& tdcOrdinal, int& totalTdcs) {
    if (!starts_with(line, "TDC ")) {
        return false;
    }
    const char* cursor = line.c_str() + 4;
    char* end = nullptr;
    const long ordinal = std::strtol(cursor, &end, 10);
    if (end == cursor) {
        return false;
    }
    cursor = end;
    while (*cursor == ' ') {
        ++cursor;
    }
    if (!starts_with(cursor, "of total ")) {
        return false;
    }
    cursor += 9;
    const long total = std::strtol(cursor, &end, 10);
    if (end == cursor) {
        return false;
    }
    tdcOrdinal = static_cast<int>(ordinal);
    totalTdcs = static_cast<int>(total);
    return true;
}

std::optional<double> parse_trigger(const std::string& line) {
    constexpr std::string_view prefix = "GlobalTriggerTime ";
    if (!starts_with(line, prefix)) {
        return std::nullopt;
    }
    const char* cursor = line.c_str() + prefix.size();
    char* end = nullptr;
    const double value = std::strtod(cursor, &end);
    if (end == cursor) {
        return std::nullopt;
    }
    return value;
}

bool parse_row(const std::string& line, int& tdcIndexOneBased, int& channel, int& isRising, double& timeNs) {
    const char* cursor = line.c_str();
    char* end = nullptr;
    const long first = std::strtol(cursor, &end, 10);
    if (end == cursor || *end != ' ') {
        return false;
    }
    cursor = end + 1;
    const long second = std::strtol(cursor, &end, 10);
    if (end == cursor || *end != ' ') {
        return false;
    }
    cursor = end + 1;
    const long third = std::strtol(cursor, &end, 10);
    if (end == cursor || *end != ' ') {
        return false;
    }
    cursor = end + 1;
    const double fourth = std::strtod(cursor, &end);
    if (end == cursor || *end != '\0') {
        return false;
    }
    tdcIndexOneBased = static_cast<int>(first);
    channel = static_cast<int>(second);
    isRising = static_cast<int>(third);
    timeNs = fourth;
    return true;
}

std::int64_t floor_div(std::int64_t numerator, std::int64_t denominator) {
    std::int64_t quotient = numerator / denominator;
    const std::int64_t remainder = numerator % denominator;
    if (remainder != 0 && numerator < 0) {
        --quotient;
    }
    return quotient;
}

std::optional<std::int64_t> select_trigger_ns(const EventState& event) {
    if (event.hasTrigger[0]) {
        return event.triggerNs[0];
    }
    if (event.hasTrigger[kTargetTdcIndex]) {
        return event.triggerNs[kTargetTdcIndex];
    }
    return std::nullopt;
}

bool is_valid_hit_time_ns(double timeNs) {
    return std::isfinite(timeNs) && std::abs(timeNs) <= kMaxAcceptedHitTimeNs;
}

void record_hit(std::unordered_map<std::int64_t, BinCounts>& bins,
                AnalysisStats& stats,
                std::int64_t binIndex,
                std::int64_t absoluteHitNs,
                const char* detector) {
    BinCounts& entry = bins[binIndex];
    if (std::string_view(detector) == "ncal1") {
        entry.ncal1 += 1;
        stats.validNcal1Hits += 1;
    } else if (std::string_view(detector) == "lstilbene") {
        entry.lStilbene += 1;
        stats.validLStilbeneHits += 1;
    } else if (std::string_view(detector) == "sstilbene") {
        entry.sStilbene += 1;
        stats.validSStilbeneHits += 1;
    }

    stats.minBinIndex = std::min(stats.minBinIndex, binIndex);
    stats.maxBinIndex = std::max(stats.maxBinIndex, binIndex);
    stats.minHitAbsoluteNs = std::min(stats.minHitAbsoluteNs, absoluteHitNs);
    stats.maxHitAbsoluteNs = std::max(stats.maxHitAbsoluteNs, absoluteHitNs);
}

void finalize_event(const Config& config,
                    AnalysisStats& stats,
                    EventState& event,
                    std::unordered_map<std::int64_t, BinCounts>& bins) {
    if (!event.hasAnyData) {
        return;
    }
    stats.finalizedEvents += 1;

    const auto triggerNs = select_trigger_ns(event);
    if (!triggerNs.has_value() || !stats.hasRunStartTrigger) {
        event = EventState{};
        return;
    }

    const std::int64_t binWidthNs = static_cast<std::int64_t>(std::llround(config.binWidthNs));
    auto handle_detector = [&](int channel, const char* detectorName) {
        const double riseTimeNs = event.times[kTargetTdcIndex][channel][1];
        const double fallTimeNs = event.times[kTargetTdcIndex][channel][0];
        if (riseTimeNs == 0.0 || fallTimeNs == 0.0) {
            return;
        }
        if (!is_valid_hit_time_ns(riseTimeNs) || !is_valid_hit_time_ns(fallTimeNs)) {
            return;
        }
        const std::int64_t absoluteHitNs = *triggerNs + static_cast<std::int64_t>(std::llround(riseTimeNs));
        const std::int64_t binIndex = floor_div(absoluteHitNs - stats.runStartTriggerNs, binWidthNs);
        record_hit(bins, stats, binIndex, absoluteHitNs, detectorName);
    };

    handle_detector(kNcal1Channel, "ncal1");
    handle_detector(kLargeStilbeneChannel, "lstilbene");
    handle_detector(kSmallStilbeneChannel, "sstilbene");

    event = EventState{};
}

void write_rates_file(const std::string& path,
                      const Config& config,
                      const AnalysisStats& stats,
                      const std::unordered_map<std::int64_t, BinCounts>& bins) {
    std::vector<std::pair<std::int64_t, BinCounts>> sortedBins;
    sortedBins.reserve(bins.size());
    for (const auto& entry : bins) {
        sortedBins.push_back(entry);
    }
    std::sort(sortedBins.begin(), sortedBins.end(), [](const auto& left, const auto& right) {
        return left.first < right.first;
    });

    std::ofstream output(path);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write rates file: " + path);
    }
    output << "# bin_width_ns=" << format_double(config.binWidthNs) << '\n';
    output << "# x_axis=time_since_run_start_s\n";
    output << "# y_axis=counts_per_bin\n";
    output << "# ncal_series=Ncal1\n";
    output << "# lstilbene_series=Lstilbene\n";
    output << "# sstilbene_series=Sstilbene\n";
    output << "# run_start_trigger_ns=" << stats.runStartTriggerNs << '\n';
    output << "# min_bin_index=" << stats.minBinIndex << '\n';
    output << "# max_bin_index=" << stats.maxBinIndex << '\n';
    output << "# columns: bin_index ncal1_count lstilbene_count sstilbene_count\n";
    for (const auto& [binIndex, counts] : sortedBins) {
        output << binIndex << ' ' << counts.ncal1 << ' ' << counts.lStilbene << ' ' << counts.sStilbene << '\n';
    }
}

void write_summary(const std::string& path,
                   const Config& config,
                   const AnalysisStats& stats,
                   std::size_t nonzeroBins,
                   double elapsedSeconds) {
    std::ofstream output(path);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write summary file: " + path);
    }
    output << "input_file=" << config.inputPath << '\n';
    output << "bin_width_ns=" << format_double(config.binWidthNs) << '\n';
    output << "trigger_reset_threshold_seconds=" << format_double(config.triggerResetThresholdSeconds) << '\n';
    output << "parsed_headers=" << stats.parsedHeaders << '\n';
    output << "parsed_rows=" << stats.parsedRows << '\n';
    output << "finalized_events=" << stats.finalizedEvents << '\n';
    output << "duplicate_restart_detected=" << (stats.duplicateRestartDetected ? "true" : "false") << '\n';
    output << "duplicate_restart_previous_trigger_seconds=" << format_double(stats.duplicateRestartPreviousTriggerSeconds) << '\n';
    output << "duplicate_restart_new_trigger_seconds=" << format_double(stats.duplicateRestartNewTriggerSeconds) << '\n';
    output << "duplicate_restart_headers_skipped=" << stats.duplicateRestartHeadersSkipped << '\n';
    output << "run_start_trigger_seconds=" << format_double(static_cast<double>(stats.runStartTriggerNs) * 1e-9) << '\n';
    output << "valid_ncal1_hits=" << stats.validNcal1Hits << '\n';
    output << "valid_lstilbene_hits=" << stats.validLStilbeneHits << '\n';
    output << "valid_sstilbene_hits=" << stats.validSStilbeneHits << '\n';
    output << "nonzero_bins=" << nonzeroBins << '\n';
    if (stats.minHitAbsoluteNs != std::numeric_limits<std::int64_t>::max() &&
        stats.maxHitAbsoluteNs != std::numeric_limits<std::int64_t>::min()) {
        output << "first_hit_time_since_run_start_seconds="
               << format_double(static_cast<double>(stats.minHitAbsoluteNs - stats.runStartTriggerNs) * 1e-9) << '\n';
        output << "last_hit_time_since_run_start_seconds="
               << format_double(static_cast<double>(stats.maxHitAbsoluteNs - stats.runStartTriggerNs) * 1e-9) << '\n';
        output << "covered_duration_seconds="
               << format_double(static_cast<double>(stats.maxHitAbsoluteNs - stats.minHitAbsoluteNs) * 1e-9) << '\n';
    }
    output << "elapsed_seconds=" << format_double(elapsedSeconds) << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Config config = parse_args(argc, argv);
        std::ifstream input(config.inputPath);
        if (!input.is_open()) {
            throw std::runtime_error("Unable to open input file: " + config.inputPath);
        }

        AnalysisStats stats;
        EventState currentEvent;
        std::unordered_map<std::int64_t, BinCounts> bins;
        bins.reserve(1 << 20);

        int currentTdcOrdinal = -1;
        bool haveLastEventStartTrigger = false;
        double lastEventStartTriggerSeconds = 0.0;
        std::string line;
        const auto startedAt = std::chrono::steady_clock::now();

        while (std::getline(input, line)) {
            int tdcOrdinal = 0;
            int totalTdcs = 0;
            if (parse_header(line, tdcOrdinal, totalTdcs)) {
                stats.parsedHeaders += 1;
                if (tdcOrdinal == 0) {
                    finalize_event(config, stats, currentEvent, bins);
                }
                currentTdcOrdinal = tdcOrdinal;
                continue;
            }

            if (const auto trigger = parse_trigger(line); trigger.has_value()) {
                if (currentTdcOrdinal >= 0 && currentTdcOrdinal < 2) {
                    const std::int64_t triggerNs = static_cast<std::int64_t>(std::llround(*trigger * 1e9));
                    currentEvent.triggerNs[static_cast<std::size_t>(currentTdcOrdinal)] = triggerNs;
                    currentEvent.hasTrigger[static_cast<std::size_t>(currentTdcOrdinal)] = true;
                    if (currentTdcOrdinal == 0) {
                        if (haveLastEventStartTrigger &&
                            *trigger + config.triggerResetThresholdSeconds < lastEventStartTriggerSeconds) {
                            stats.duplicateRestartDetected = true;
                            stats.duplicateRestartPreviousTriggerSeconds = lastEventStartTriggerSeconds;
                            stats.duplicateRestartNewTriggerSeconds = *trigger;
                            stats.duplicateRestartHeadersSkipped = stats.parsedHeaders;
                            break;
                        }
                        if (!stats.hasRunStartTrigger) {
                            stats.hasRunStartTrigger = true;
                            stats.runStartTriggerNs = triggerNs;
                        }
                        haveLastEventStartTrigger = true;
                        lastEventStartTriggerSeconds = *trigger;
                    }
                }
                continue;
            }

            int tdcIndexOneBased = 0;
            int channel = 0;
            int isRising = 0;
            double timeNs = 0.0;
            if (!parse_row(line, tdcIndexOneBased, channel, isRising, timeNs)) {
                continue;
            }
            stats.parsedRows += 1;
            if (currentTdcOrdinal < 0 || currentTdcOrdinal >= 2) {
                continue;
            }
            if (channel < 0 || channel >= 32) {
                continue;
            }
            if (isRising != 0 && isRising != 1) {
                continue;
            }

            currentEvent.times[static_cast<std::size_t>(currentTdcOrdinal)][static_cast<std::size_t>(channel)][static_cast<std::size_t>(isRising)] = timeNs;
            currentEvent.hasAnyData = true;
        }

        finalize_event(config, stats, currentEvent, bins);

        const auto endedAt = std::chrono::steady_clock::now();
        const double elapsedSeconds = std::chrono::duration_cast<std::chrono::duration<double>>(endedAt - startedAt).count();

        write_rates_file(config.outputPrefix + "_rates.txt", config, stats, bins);
        write_summary(config.outputPrefix + "_summary.txt", config, stats, bins.size(), elapsedSeconds);

        std::cout << "Input file: " << config.inputPath << '\n';
        std::cout << "Output prefix: " << config.outputPrefix << '\n';
        std::cout << "Finalized events: " << stats.finalizedEvents << '\n';
        std::cout << "Duplicate restart detected: " << (stats.duplicateRestartDetected ? "yes" : "no") << '\n';
        std::cout << "Valid Ncal1 hits: " << stats.validNcal1Hits << '\n';
        std::cout << "Valid Lstilbene hits: " << stats.validLStilbeneHits << '\n';
        std::cout << "Valid Sstilbene hits: " << stats.validSStilbeneHits << '\n';
        std::cout << "Nonzero 1 us bins: " << bins.size() << '\n';
        std::cout << "Elapsed seconds: " << format_double(elapsedSeconds) << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}