#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace {

struct Config {
    std::string inputPath;
    std::string outputPrefix;
    double triggerResetThresholdSeconds = 1.0;
};

struct Histogram1D {
    std::string objectName;
    std::string title;
    std::string xTitle;
    int bins = 0;
    double min = 0.0;
    double max = 0.0;
    std::vector<std::uint64_t> storage;
    double entries = 0.0;
    double sumw = 0.0;
    double sumw2 = 0.0;
    double sumwx = 0.0;
    double sumwx2 = 0.0;
};

struct AnalysisStats {
    std::uint64_t parsedHeaders = 0;
    std::uint64_t parsedRows = 0;
    std::uint64_t validNcal1Rises = 0;
    std::uint64_t validLStilbeneRises = 0;
    std::uint64_t validSStilbeneRises = 0;
    bool duplicateRestartDetected = false;
    double duplicateRestartPreviousTriggerSeconds = 0.0;
    double duplicateRestartNewTriggerSeconds = 0.0;
    std::uint64_t duplicateRestartHeadersSkipped = 0;
};

constexpr int kTargetTdcOrdinal = 1;
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
        } else if (arg == "--trigger-reset-threshold-seconds" && index + 1 < argc) {
            config.triggerResetThresholdSeconds = std::stod(argv[++index]);
        } else if (arg == "--help") {
            std::cout
                << "Usage: dogma_trgref_artifact_hunting --input <file> --output-prefix <prefix>\n"
                << "  [--trigger-reset-threshold-seconds 1.0]\n";
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

Histogram1D make_histogram_1d(const std::string& objectName,
                              const std::string& title,
                              int bins,
                              double min,
                              double max,
                              const std::string& xTitle) {
    Histogram1D histogram;
    histogram.objectName = objectName;
    histogram.title = title;
    histogram.xTitle = xTitle;
    histogram.bins = bins;
    histogram.min = min;
    histogram.max = max;
    histogram.storage.assign(static_cast<std::size_t>(bins) + 2, 0);
    return histogram;
}

void fill_histogram(Histogram1D& histogram, double value) {
    histogram.entries += 1.0;
    histogram.sumw += 1.0;
    histogram.sumw2 += 1.0;
    histogram.sumwx += value;
    histogram.sumwx2 += value * value;

    int binIndex = 0;
    if (value >= histogram.max) {
        binIndex = histogram.bins + 1;
    } else if (value >= histogram.min) {
        const double fraction = (value - histogram.min) / (histogram.max - histogram.min);
        binIndex = 1 + std::min(histogram.bins - 1, static_cast<int>(std::floor(fraction * histogram.bins)));
    }
    histogram.storage[static_cast<std::size_t>(binIndex)] += 1;
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

bool is_valid_hit_time_ns(double timeNs) {
    return std::isfinite(timeNs) && std::abs(timeNs) <= kMaxAcceptedHitTimeNs;
}

void write_histogram_1d(const std::string& path, const Histogram1D& histogram) {
    std::ofstream output(path);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write histogram file: " + path);
    }
    output << "# object_name=" << histogram.objectName << '\n';
    output << "# title=" << histogram.title << '\n';
    output << "# x_title=" << histogram.xTitle << '\n';
    output << "# bins=" << histogram.bins << '\n';
    output << "# min=" << format_double(histogram.min) << '\n';
    output << "# max=" << format_double(histogram.max) << '\n';
    output << "# entries=" << format_double(histogram.entries) << '\n';
    output << "# sumw=" << format_double(histogram.sumw) << '\n';
    output << "# sumw2=" << format_double(histogram.sumw2) << '\n';
    output << "# sumwx=" << format_double(histogram.sumwx) << '\n';
    output << "# sumwx2=" << format_double(histogram.sumwx2) << '\n';
    output << "# columns: bin_index count\n";
    for (std::size_t index = 0; index < histogram.storage.size(); ++index) {
        const auto count = histogram.storage[index];
        if (count == 0) {
            continue;
        }
        output << index << ' ' << count << '\n';
    }
}

void write_summary(const std::string& path,
                   const Config& config,
                   const AnalysisStats& stats,
                   const Histogram1D& hNcal1,
                   const Histogram1D& hLStilbene,
                   const Histogram1D& hSStilbene,
                   double elapsedSeconds) {
    std::ofstream output(path);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write summary file: " + path);
    }
    output << "input_file=" << config.inputPath << '\n';
    output << "trigger_reset_threshold_seconds=" << format_double(config.triggerResetThresholdSeconds) << '\n';
    output << "parsed_headers=" << stats.parsedHeaders << '\n';
    output << "parsed_rows=" << stats.parsedRows << '\n';
    output << "duplicate_restart_detected=" << (stats.duplicateRestartDetected ? "true" : "false") << '\n';
    output << "duplicate_restart_previous_trigger_seconds=" << format_double(stats.duplicateRestartPreviousTriggerSeconds) << '\n';
    output << "duplicate_restart_new_trigger_seconds=" << format_double(stats.duplicateRestartNewTriggerSeconds) << '\n';
    output << "duplicate_restart_headers_skipped=" << stats.duplicateRestartHeadersSkipped << '\n';
    output << "valid_ncal1_rises=" << stats.validNcal1Rises << '\n';
    output << "valid_lstilbene_rises=" << stats.validLStilbeneRises << '\n';
    output << "valid_sstilbene_rises=" << stats.validSStilbeneRises << '\n';
    output << "hNcal1_entries=" << format_double(hNcal1.entries) << '\n';
    output << "hLStilbene_entries=" << format_double(hLStilbene.entries) << '\n';
    output << "hSStilbene_entries=" << format_double(hSStilbene.entries) << '\n';
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

        Histogram1D hNcal1 = make_histogram_1d(
            "custom_RawRiseCountsNcal1Folded",
            "Ncal1 raw rising-edge counts vs trigger-relative time",
            12000,
            -6000.0,
            6000.0,
            "T (ns)"
        );
        Histogram1D hLStilbene = make_histogram_1d(
            "custom_RawRiseCountsLstilbeneFolded",
            "Lstilbene raw rising-edge counts vs trigger-relative time",
            12000,
            -6000.0,
            6000.0,
            "T (ns)"
        );
        Histogram1D hSStilbene = make_histogram_1d(
            "custom_RawRiseCountsSstilbeneFolded",
            "Sstilbene raw rising-edge counts vs trigger-relative time",
            12000,
            -6000.0,
            6000.0,
            "T (ns)"
        );

        AnalysisStats stats;
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
                currentTdcOrdinal = tdcOrdinal;
                continue;
            }

            if (const auto trigger = parse_trigger(line); trigger.has_value()) {
                if (currentTdcOrdinal == 0) {
                    if (haveLastEventStartTrigger &&
                        *trigger + config.triggerResetThresholdSeconds < lastEventStartTriggerSeconds) {
                        stats.duplicateRestartDetected = true;
                        stats.duplicateRestartPreviousTriggerSeconds = lastEventStartTriggerSeconds;
                        stats.duplicateRestartNewTriggerSeconds = *trigger;
                        stats.duplicateRestartHeadersSkipped = stats.parsedHeaders;
                        break;
                    }
                    haveLastEventStartTrigger = true;
                    lastEventStartTriggerSeconds = *trigger;
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
            if (currentTdcOrdinal != kTargetTdcOrdinal) {
                continue;
            }
            if (isRising != 1) {
                continue;
            }
            if (!is_valid_hit_time_ns(timeNs)) {
                continue;
            }

            if (channel == kNcal1Channel) {
                fill_histogram(hNcal1, timeNs);
                stats.validNcal1Rises += 1;
            } else if (channel == kLargeStilbeneChannel) {
                fill_histogram(hLStilbene, timeNs);
                stats.validLStilbeneRises += 1;
            } else if (channel == kSmallStilbeneChannel) {
                fill_histogram(hSStilbene, timeNs);
                stats.validSStilbeneRises += 1;
            }
        }

        const auto endedAt = std::chrono::steady_clock::now();
        const double elapsedSeconds = std::chrono::duration_cast<std::chrono::duration<double>>(endedAt - startedAt).count();

        write_histogram_1d(config.outputPrefix + "_custom_RawRiseCountsNcal1Folded_hist.txt", hNcal1);
        write_histogram_1d(config.outputPrefix + "_custom_RawRiseCountsLstilbeneFolded_hist.txt", hLStilbene);
        write_histogram_1d(config.outputPrefix + "_custom_RawRiseCountsSstilbeneFolded_hist.txt", hSStilbene);
        write_summary(config.outputPrefix + "_summary.txt", config, stats, hNcal1, hLStilbene, hSStilbene, elapsedSeconds);

        std::cout << "Input file: " << config.inputPath << '\n';
        std::cout << "Output prefix: " << config.outputPrefix << '\n';
        std::cout << "Duplicate restart detected: " << (stats.duplicateRestartDetected ? "yes" : "no") << '\n';
        std::cout << "Ncal1 raw rises: " << stats.validNcal1Rises << '\n';
        std::cout << "Lstilbene raw rises: " << stats.validLStilbeneRises << '\n';
        std::cout << "Sstilbene raw rises: " << stats.validSStilbeneRises << '\n';
        std::cout << "Elapsed seconds: " << format_double(elapsedSeconds) << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}