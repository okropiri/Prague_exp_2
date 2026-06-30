#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
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
    std::uint64_t validRiseEdges = 0;
    std::uint64_t validFallEdges = 0;
    std::uint64_t windowsWithRise = 0;
    std::uint64_t windowsWithTwoPlusRises = 0;
    std::uint64_t windowsWithThreePlusRises = 0;
    std::uint64_t totalRiseToRisePairs = 0;
    std::uint64_t riseSpacingLt1Ns = 0;
    std::uint64_t riseSpacingLt2Ns = 0;
    std::uint64_t riseSpacingLt5Ns = 0;
    std::uint64_t riseSpacingLt10Ns = 0;
    std::uint64_t riseSpacingLt20Ns = 0;
    bool duplicateRestartDetected = false;
    double duplicateRestartPreviousTriggerSeconds = 0.0;
    double duplicateRestartNewTriggerSeconds = 0.0;
    std::uint64_t duplicateRestartHeadersSkipped = 0;
};

constexpr int kTargetTdcOrdinal = 1;
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
                << "Usage: dogma_lstilbene_edge_checks --input <file> --output-prefix <prefix>\n"
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

double mean_bin_content_in_range(const Histogram1D& histogram, double minValue, double maxValue) {
    const double binWidth = (histogram.max - histogram.min) / static_cast<double>(histogram.bins);
    double sum = 0.0;
    std::uint64_t count = 0;
    for (int bin = 1; bin <= histogram.bins; ++bin) {
        const double center = histogram.min + (static_cast<double>(bin) - 0.5) * binWidth;
        if (center < minValue || center >= maxValue) {
            continue;
        }
        sum += static_cast<double>(histogram.storage[static_cast<std::size_t>(bin)]);
        count += 1;
    }
    return count == 0 ? 0.0 : sum / static_cast<double>(count);
}

void finalize_lstilbene_window(std::vector<double>& riseTimes,
                               Histogram1D& hRiseSpacing,
                               AnalysisStats& stats) {
    if (riseTimes.empty()) {
        return;
    }

    stats.windowsWithRise += 1;
    if (riseTimes.size() >= 2) {
        stats.windowsWithTwoPlusRises += 1;
    }
    if (riseTimes.size() >= 3) {
        stats.windowsWithThreePlusRises += 1;
    }

    std::sort(riseTimes.begin(), riseTimes.end());
    for (std::size_t index = 1; index < riseTimes.size(); ++index) {
        const double dt = riseTimes[index] - riseTimes[index - 1];
        if (!std::isfinite(dt) || dt < 0.0) {
            continue;
        }
        fill_histogram(hRiseSpacing, dt);
        stats.totalRiseToRisePairs += 1;
        if (dt < 1.0) {
            stats.riseSpacingLt1Ns += 1;
        }
        if (dt < 2.0) {
            stats.riseSpacingLt2Ns += 1;
        }
        if (dt < 5.0) {
            stats.riseSpacingLt5Ns += 1;
        }
        if (dt < 10.0) {
            stats.riseSpacingLt10Ns += 1;
        }
        if (dt < 20.0) {
            stats.riseSpacingLt20Ns += 1;
        }
    }
    riseTimes.clear();
}

void write_summary(const std::string& path,
                   const Config& config,
                   const AnalysisStats& stats,
                   const Histogram1D& hRiseFolded,
                   const Histogram1D& hFallFolded,
                   const Histogram1D& hRiseSpacing,
                   double elapsedSeconds) {
    const double riseLeftMean = mean_bin_content_in_range(hRiseFolded, -5000.0, -2500.0);
    const double riseMidMean = mean_bin_content_in_range(hRiseFolded, -2500.0, 2500.0);
    const double riseRightMean = mean_bin_content_in_range(hRiseFolded, 2500.0, 5000.0);
    const double fallLeftMean = mean_bin_content_in_range(hFallFolded, -5000.0, -2500.0);
    const double fallMidMean = mean_bin_content_in_range(hFallFolded, -2500.0, 2500.0);
    const double fallRightMean = mean_bin_content_in_range(hFallFolded, 2500.0, 5000.0);

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
    output << "valid_lstilbene_rises=" << stats.validRiseEdges << '\n';
    output << "valid_lstilbene_falls=" << stats.validFallEdges << '\n';
    output << "lstilbene_windows_with_rise=" << stats.windowsWithRise << '\n';
    output << "lstilbene_windows_with_two_plus_rises=" << stats.windowsWithTwoPlusRises << '\n';
    output << "lstilbene_windows_with_three_plus_rises=" << stats.windowsWithThreePlusRises << '\n';
    output << "rise_to_rise_pairs=" << stats.totalRiseToRisePairs << '\n';
    output << "rise_spacing_lt_1ns=" << stats.riseSpacingLt1Ns << '\n';
    output << "rise_spacing_lt_2ns=" << stats.riseSpacingLt2Ns << '\n';
    output << "rise_spacing_lt_5ns=" << stats.riseSpacingLt5Ns << '\n';
    output << "rise_spacing_lt_10ns=" << stats.riseSpacingLt10Ns << '\n';
    output << "rise_spacing_lt_20ns=" << stats.riseSpacingLt20Ns << '\n';
    output << "rise_folded_left_mean=" << format_double(riseLeftMean) << '\n';
    output << "rise_folded_mid_mean=" << format_double(riseMidMean) << '\n';
    output << "rise_folded_right_mean=" << format_double(riseRightMean) << '\n';
    output << "rise_folded_right_to_left_ratio="
           << format_double(riseLeftMean == 0.0 ? 0.0 : riseRightMean / riseLeftMean) << '\n';
    output << "fall_folded_left_mean=" << format_double(fallLeftMean) << '\n';
    output << "fall_folded_mid_mean=" << format_double(fallMidMean) << '\n';
    output << "fall_folded_right_mean=" << format_double(fallRightMean) << '\n';
    output << "fall_folded_right_to_left_ratio="
           << format_double(fallLeftMean == 0.0 ? 0.0 : fallRightMean / fallLeftMean) << '\n';
    output << "rise_spacing_entries=" << format_double(hRiseSpacing.entries) << '\n';
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

        Histogram1D hRiseFolded = make_histogram_1d(
            "custom_RawRiseCountsLstilbeneFolded",
            "Lstilbene raw rising-edge counts vs trigger-relative time",
            12000,
            -6000.0,
            6000.0,
            "T (ns)"
        );
        Histogram1D hFallFolded = make_histogram_1d(
            "custom_RawFallCountsLstilbeneFolded",
            "Lstilbene raw falling-edge counts vs trigger-relative time",
            12000,
            -6000.0,
            6000.0,
            "T (ns)"
        );
        Histogram1D hRiseSpacing = make_histogram_1d(
            "custom_RiseToRiseSpacingLstilbene",
            "Lstilbene rise-to-rise spacing within trigger windows",
            2000,
            0.0,
            1000.0,
            "Delta T between successive rises (ns)"
        );

        AnalysisStats stats;
        int currentTdcOrdinal = -1;
        bool haveLastEventStartTrigger = false;
        double lastEventStartTriggerSeconds = 0.0;
        std::vector<double> currentWindowRiseTimes;
        std::string line;
        const auto startedAt = std::chrono::steady_clock::now();

        while (std::getline(input, line)) {
            int tdcOrdinal = 0;
            int totalTdcs = 0;
            if (parse_header(line, tdcOrdinal, totalTdcs)) {
                if (currentTdcOrdinal == kTargetTdcOrdinal) {
                    finalize_lstilbene_window(currentWindowRiseTimes, hRiseSpacing, stats);
                }
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
            if (channel != kLargeStilbeneChannel) {
                continue;
            }
            if (!is_valid_hit_time_ns(timeNs)) {
                continue;
            }

            if (isRising == 1) {
                fill_histogram(hRiseFolded, timeNs);
                currentWindowRiseTimes.push_back(timeNs);
                stats.validRiseEdges += 1;
            } else if (isRising == 0) {
                fill_histogram(hFallFolded, timeNs);
                stats.validFallEdges += 1;
            }
        }

        if (currentTdcOrdinal == kTargetTdcOrdinal) {
            finalize_lstilbene_window(currentWindowRiseTimes, hRiseSpacing, stats);
        }

        const auto endedAt = std::chrono::steady_clock::now();
        const double elapsedSeconds = std::chrono::duration_cast<std::chrono::duration<double>>(endedAt - startedAt).count();

        write_histogram_1d(config.outputPrefix + "_custom_RawRiseCountsLstilbeneFolded_hist.txt", hRiseFolded);
        write_histogram_1d(config.outputPrefix + "_custom_RawFallCountsLstilbeneFolded_hist.txt", hFallFolded);
        write_histogram_1d(config.outputPrefix + "_custom_RiseToRiseSpacingLstilbene_hist.txt", hRiseSpacing);
        write_summary(config.outputPrefix + "_summary.txt", config, stats, hRiseFolded, hFallFolded, hRiseSpacing, elapsedSeconds);

        std::cout << "Lstilbene rise edges: " << stats.validRiseEdges << '\n';
        std::cout << "Lstilbene fall edges: " << stats.validFallEdges << '\n';
        std::cout << "Rise-to-rise pairs: " << stats.totalRiseToRisePairs << '\n';
        std::cout << "Duplicate restart detected: " << (stats.duplicateRestartDetected ? "true" : "false") << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}