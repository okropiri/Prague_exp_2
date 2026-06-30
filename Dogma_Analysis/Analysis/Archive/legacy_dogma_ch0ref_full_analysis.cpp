#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <deque>
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
#include <vector>

namespace {

struct Config {
    std::string inputPath;
    std::string outputPrefix;
    double foldedMinNs = -6000.0;
    double foldedMaxNs = 6000.0;
    int foldedBins = 12000;
    double totMinNs = 0.0;
    double totMaxNs = 128.0;
    int totBins = 512;
    double rateBinWidthNs = 1000.0;
    std::uint64_t asymmetryBlockSizeWindows = 1000;
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

struct Histogram2D {
    std::string objectName;
    std::string title;
    std::string xTitle;
    std::string yTitle;
    int xBins = 0;
    double xMin = 0.0;
    double xMax = 0.0;
    int yBins = 0;
    double yMin = 0.0;
    double yMax = 0.0;
    std::vector<std::uint64_t> storage;
    double entries = 0.0;
    double sumw = 0.0;
    double sumw2 = 0.0;
    double sumwx = 0.0;
    double sumwx2 = 0.0;
    double sumwy = 0.0;
    double sumwy2 = 0.0;
    double sumwxy = 0.0;
};

struct RateBinCounts {
    std::uint32_t ch0 = 0;
    std::uint32_t ncal1 = 0;
    std::uint32_t lStilbene = 0;
    std::uint32_t sStilbene = 0;
};

struct WindowDetectorCounts {
    std::uint32_t total = 0;
    std::uint32_t after = 0;
};

struct AggregateCounts {
    std::uint64_t total = 0;
    std::uint64_t after = 0;
    std::uint64_t nonzeroWindows = 0;
};

struct AsymmetryBlock {
    std::uint64_t windowCount = 0;
    std::uint64_t validCh0Windows = 0;
    AggregateCounts ncal1;
    AggregateCounts lStilbene;
    AggregateCounts sStilbene;
};

struct PulsePair {
    double riseNs = 0.0;
    double totNs = 0.0;
};

struct WindowState {
    bool hasTrigger = false;
    std::int64_t triggerNs = 0;
    std::vector<double> ch0Rises;
    std::vector<double> ncal1Rises;
    std::vector<double> lStilbeneRises;
    std::vector<double> sStilbeneRises;
    std::deque<double> pendingNcal1Rises;
    std::deque<double> pendingLStilbeneRises;
    std::deque<double> pendingSStilbeneRises;
    std::vector<PulsePair> ncal1Pairs;
    std::vector<PulsePair> lStilbenePairs;
    std::vector<PulsePair> sStilbenePairs;
};

struct AnalysisStats {
    std::uint64_t parsedHeaders = 0;
    std::uint64_t parsedRows = 0;
    std::uint64_t finalizedWindows = 0;
    std::uint64_t validCh0Rises = 0;
    std::uint64_t validNcal1Rises = 0;
    std::uint64_t validLStilbeneRises = 0;
    std::uint64_t validSStilbeneRises = 0;
    std::uint64_t validNcal1Pairs = 0;
    std::uint64_t validLStilbenePairs = 0;
    std::uint64_t validSStilbenePairs = 0;
    std::uint64_t windowsWithValidCh0 = 0;
    std::uint64_t windowsWithoutCh0 = 0;
    std::uint64_t windowsWithMultipleCh0 = 0;
    bool duplicateRestartDetected = false;
    double duplicateRestartPreviousTriggerSeconds = 0.0;
    double duplicateRestartNewTriggerSeconds = 0.0;
    std::uint64_t duplicateRestartHeadersSkipped = 0;
    bool hasRunStartTrigger = false;
    std::int64_t runStartTriggerNs = 0;
    std::int64_t minHitAbsoluteNs = std::numeric_limits<std::int64_t>::max();
    std::int64_t maxHitAbsoluteNs = std::numeric_limits<std::int64_t>::min();
    std::int64_t minRateBinIndex = std::numeric_limits<std::int64_t>::max();
    std::int64_t maxRateBinIndex = std::numeric_limits<std::int64_t>::min();
    AggregateCounts ncal1Asymmetry;
    AggregateCounts lStilbeneAsymmetry;
    AggregateCounts sStilbeneAsymmetry;
};

constexpr int kTriggerOrdinal = 0;
constexpr int kSignalOrdinal = 1;
constexpr int kCh0Channel = 0;
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
        } else if (arg == "--folded-min-ns" && index + 1 < argc) {
            config.foldedMinNs = std::stod(argv[++index]);
        } else if (arg == "--folded-max-ns" && index + 1 < argc) {
            config.foldedMaxNs = std::stod(argv[++index]);
        } else if (arg == "--folded-bins" && index + 1 < argc) {
            config.foldedBins = std::stoi(argv[++index]);
        } else if (arg == "--tot-min-ns" && index + 1 < argc) {
            config.totMinNs = std::stod(argv[++index]);
        } else if (arg == "--tot-max-ns" && index + 1 < argc) {
            config.totMaxNs = std::stod(argv[++index]);
        } else if (arg == "--tot-bins" && index + 1 < argc) {
            config.totBins = std::stoi(argv[++index]);
        } else if (arg == "--rate-bin-width-ns" && index + 1 < argc) {
            config.rateBinWidthNs = std::stod(argv[++index]);
        } else if (arg == "--asymmetry-block-size-windows" && index + 1 < argc) {
            config.asymmetryBlockSizeWindows = static_cast<std::uint64_t>(std::stoull(argv[++index]));
        } else if (arg == "--trigger-reset-threshold-seconds" && index + 1 < argc) {
            config.triggerResetThresholdSeconds = std::stod(argv[++index]);
        } else if (arg == "--help") {
            std::cout
                << "Usage: dogma_ch0ref_full_analysis --input <file> --output-prefix <prefix>\n"
                << "  [--folded-min-ns -6000] [--folded-max-ns 6000] [--folded-bins 12000]\n"
                << "  [--tot-min-ns 0] [--tot-max-ns 128] [--tot-bins 512]\n"
                << "  [--rate-bin-width-ns 1000] [--asymmetry-block-size-windows 1000]\n"
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
    if (config.foldedBins <= 0 || config.totBins <= 0) {
        throw std::runtime_error("Histogram bin counts must be positive");
    }
    if (!(config.foldedMinNs < config.foldedMaxNs)) {
        throw std::runtime_error("folded-min-ns must be smaller than folded-max-ns");
    }
    if (!(config.totMinNs < config.totMaxNs)) {
        throw std::runtime_error("tot-min-ns must be smaller than tot-max-ns");
    }
    if (config.rateBinWidthNs <= 0.0) {
        throw std::runtime_error("rate-bin-width-ns must be positive");
    }
    if (config.asymmetryBlockSizeWindows == 0) {
        throw std::runtime_error("asymmetry-block-size-windows must be positive");
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

Histogram2D make_histogram_2d(const std::string& objectName,
                              const std::string& title,
                              int xBins,
                              double xMin,
                              double xMax,
                              const std::string& xTitle,
                              int yBins,
                              double yMin,
                              double yMax,
                              const std::string& yTitle) {
    Histogram2D histogram;
    histogram.objectName = objectName;
    histogram.title = title;
    histogram.xTitle = xTitle;
    histogram.yTitle = yTitle;
    histogram.xBins = xBins;
    histogram.xMin = xMin;
    histogram.xMax = xMax;
    histogram.yBins = yBins;
    histogram.yMin = yMin;
    histogram.yMax = yMax;
    histogram.storage.assign(static_cast<std::size_t>(xBins + 2) * static_cast<std::size_t>(yBins + 2), 0);
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

void fill_histogram(Histogram2D& histogram, double xValue, double yValue) {
    histogram.entries += 1.0;
    histogram.sumw += 1.0;
    histogram.sumw2 += 1.0;
    histogram.sumwx += xValue;
    histogram.sumwx2 += xValue * xValue;
    histogram.sumwy += yValue;
    histogram.sumwy2 += yValue * yValue;
    histogram.sumwxy += xValue * yValue;

    int xBinIndex = 0;
    if (xValue >= histogram.xMax) {
        xBinIndex = histogram.xBins + 1;
    } else if (xValue >= histogram.xMin) {
        const double fraction = (xValue - histogram.xMin) / (histogram.xMax - histogram.xMin);
        xBinIndex = 1 + std::min(histogram.xBins - 1, static_cast<int>(std::floor(fraction * histogram.xBins)));
    }

    int yBinIndex = 0;
    if (yValue >= histogram.yMax) {
        yBinIndex = histogram.yBins + 1;
    } else if (yValue >= histogram.yMin) {
        const double fraction = (yValue - histogram.yMin) / (histogram.yMax - histogram.yMin);
        yBinIndex = 1 + std::min(histogram.yBins - 1, static_cast<int>(std::floor(fraction * histogram.yBins)));
    }

    const std::size_t flatIndex = static_cast<std::size_t>(yBinIndex) * static_cast<std::size_t>(histogram.xBins + 2)
                                + static_cast<std::size_t>(xBinIndex);
    histogram.storage[flatIndex] += 1;
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

std::int64_t floor_div(std::int64_t numerator, std::int64_t denominator) {
    std::int64_t quotient = numerator / denominator;
    const std::int64_t remainder = numerator % denominator;
    if (remainder != 0 && numerator < 0) {
        --quotient;
    }
    return quotient;
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

void write_histogram_2d(const std::string& path, const Histogram2D& histogram) {
    std::ofstream output(path);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write histogram file: " + path);
    }
    output << "# object_name=" << histogram.objectName << '\n';
    output << "# title=" << histogram.title << '\n';
    output << "# x_title=" << histogram.xTitle << '\n';
    output << "# y_title=" << histogram.yTitle << '\n';
    output << "# x_bins=" << histogram.xBins << '\n';
    output << "# x_min=" << format_double(histogram.xMin) << '\n';
    output << "# x_max=" << format_double(histogram.xMax) << '\n';
    output << "# y_bins=" << histogram.yBins << '\n';
    output << "# y_min=" << format_double(histogram.yMin) << '\n';
    output << "# y_max=" << format_double(histogram.yMax) << '\n';
    output << "# entries=" << format_double(histogram.entries) << '\n';
    output << "# sumw=" << format_double(histogram.sumw) << '\n';
    output << "# sumw2=" << format_double(histogram.sumw2) << '\n';
    output << "# sumwx=" << format_double(histogram.sumwx) << '\n';
    output << "# sumwx2=" << format_double(histogram.sumwx2) << '\n';
    output << "# sumwy=" << format_double(histogram.sumwy) << '\n';
    output << "# sumwy2=" << format_double(histogram.sumwy2) << '\n';
    output << "# sumwxy=" << format_double(histogram.sumwxy) << '\n';
    output << "# columns: x_bin_index y_bin_index count\n";
    for (int yIndex = 0; yIndex < histogram.yBins + 2; ++yIndex) {
        for (int xIndex = 0; xIndex < histogram.xBins + 2; ++xIndex) {
            const std::size_t flatIndex = static_cast<std::size_t>(yIndex) * static_cast<std::size_t>(histogram.xBins + 2)
                                        + static_cast<std::size_t>(xIndex);
            const auto count = histogram.storage[flatIndex];
            if (count == 0) {
                continue;
            }
            output << xIndex << ' ' << yIndex << ' ' << count << '\n';
        }
    }
}

void accumulate_counts(AggregateCounts& aggregate, const WindowDetectorCounts& counts) {
    aggregate.total += counts.total;
    aggregate.after += counts.after;
    if (counts.total > 0) {
        aggregate.nonzeroWindows += 1;
    }
}

AsymmetryBlock& ensure_block(std::vector<AsymmetryBlock>& blocks,
                             std::uint64_t windowIndex,
                             std::uint64_t blockSizeWindows) {
    const std::size_t blockIndex = static_cast<std::size_t>(windowIndex / blockSizeWindows);
    if (blocks.size() <= blockIndex) {
        blocks.resize(blockIndex + 1);
    }
    return blocks[blockIndex];
}

void write_asymmetry_blocks_file(const std::string& path,
                                 const Config& config,
                                 const std::vector<AsymmetryBlock>& blocks) {
    std::ofstream output(path);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write asymmetry blocks file: " + path);
    }
    output << "# block_size_windows=" << config.asymmetryBlockSizeWindows << '\n';
    output << "# folded_min_ns=" << format_double(config.foldedMinNs) << '\n';
    output << "# folded_max_ns=" << format_double(config.foldedMaxNs) << '\n';
    output << "# columns: block_index window_start window_end_exclusive windows_in_block valid_ch0_windows"
           << " ncal1_after ncal1_total lstilbene_after lstilbene_total sstilbene_after sstilbene_total\n";
    for (std::size_t blockIndex = 0; blockIndex < blocks.size(); ++blockIndex) {
        const AsymmetryBlock& block = blocks[blockIndex];
        const std::uint64_t windowStart = static_cast<std::uint64_t>(blockIndex) * config.asymmetryBlockSizeWindows;
        const std::uint64_t windowEndExclusive = windowStart + block.windowCount;
        output << blockIndex << ' '
               << windowStart << ' '
               << windowEndExclusive << ' '
               << block.windowCount << ' '
               << block.validCh0Windows << ' '
               << block.ncal1.after << ' ' << block.ncal1.total << ' '
               << block.lStilbene.after << ' ' << block.lStilbene.total << ' '
               << block.sStilbene.after << ' ' << block.sStilbene.total << '\n';
    }
}

void write_rates_file(const std::string& path,
                      const Config& config,
                      const AnalysisStats& stats,
                      const std::unordered_map<std::int64_t, RateBinCounts>& bins) {
    std::vector<std::pair<std::int64_t, RateBinCounts>> sortedBins;
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
    output << "# rate_bin_width_ns=" << format_double(config.rateBinWidthNs) << '\n';
    output << "# run_start_trigger_ns=" << stats.runStartTriggerNs << '\n';
    output << "# min_bin_index=" << stats.minRateBinIndex << '\n';
    output << "# max_bin_index=" << stats.maxRateBinIndex << '\n';
    output << "# columns: bin_index ch0_count ncal1_count lstilbene_count sstilbene_count\n";
    for (const auto& [binIndex, counts] : sortedBins) {
        output << binIndex << ' '
               << counts.ch0 << ' '
               << counts.ncal1 << ' '
               << counts.lStilbene << ' '
               << counts.sStilbene << '\n';
    }
}

void write_summary(const std::string& path,
                   const Config& config,
                   const AnalysisStats& stats,
                   const std::vector<AsymmetryBlock>& blocks,
                   double elapsedSeconds) {
    std::ofstream output(path);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write summary file: " + path);
    }

    const auto rateSeconds = stats.minHitAbsoluteNs <= stats.maxHitAbsoluteNs
                           ? static_cast<double>(stats.maxHitAbsoluteNs - stats.runStartTriggerNs) * 1.0e-9
                           : 0.0;
    const auto fraction = [](const AggregateCounts& counts) {
        return counts.total == 0 ? 0.0 : static_cast<double>(counts.after) / static_cast<double>(counts.total);
    };

    output << "input_file=" << config.inputPath << '\n';
    output << "folded_min_ns=" << format_double(config.foldedMinNs) << '\n';
    output << "folded_max_ns=" << format_double(config.foldedMaxNs) << '\n';
    output << "folded_bins=" << config.foldedBins << '\n';
    output << "tot_min_ns=" << format_double(config.totMinNs) << '\n';
    output << "tot_max_ns=" << format_double(config.totMaxNs) << '\n';
    output << "tot_bins=" << config.totBins << '\n';
    output << "rate_bin_width_ns=" << format_double(config.rateBinWidthNs) << '\n';
    output << "asymmetry_block_size_windows=" << config.asymmetryBlockSizeWindows << '\n';
    output << "trigger_reset_threshold_seconds=" << format_double(config.triggerResetThresholdSeconds) << '\n';
    output << "parsed_headers=" << stats.parsedHeaders << '\n';
    output << "parsed_rows=" << stats.parsedRows << '\n';
    output << "finalized_windows=" << stats.finalizedWindows << '\n';
    output << "asymmetry_block_count=" << blocks.size() << '\n';
    output << "duplicate_restart_detected=" << (stats.duplicateRestartDetected ? "true" : "false") << '\n';
    output << "duplicate_restart_previous_trigger_seconds=" << format_double(stats.duplicateRestartPreviousTriggerSeconds) << '\n';
    output << "duplicate_restart_new_trigger_seconds=" << format_double(stats.duplicateRestartNewTriggerSeconds) << '\n';
    output << "duplicate_restart_headers_skipped=" << stats.duplicateRestartHeadersSkipped << '\n';
    output << "run_start_trigger_seconds=" << format_double(static_cast<double>(stats.runStartTriggerNs) * 1.0e-9) << '\n';
    output << "last_hit_time_since_run_start_seconds=" << format_double(rateSeconds) << '\n';
    output << "valid_ch0_rises=" << stats.validCh0Rises << '\n';
    output << "valid_ncal1_rises=" << stats.validNcal1Rises << '\n';
    output << "valid_lstilbene_rises=" << stats.validLStilbeneRises << '\n';
    output << "valid_sstilbene_rises=" << stats.validSStilbeneRises << '\n';
    output << "valid_ncal1_tot_pairs=" << stats.validNcal1Pairs << '\n';
    output << "valid_lstilbene_tot_pairs=" << stats.validLStilbenePairs << '\n';
    output << "valid_sstilbene_tot_pairs=" << stats.validSStilbenePairs << '\n';
    output << "windows_with_valid_ch0=" << stats.windowsWithValidCh0 << '\n';
    output << "windows_without_ch0=" << stats.windowsWithoutCh0 << '\n';
    output << "windows_with_multiple_ch0=" << stats.windowsWithMultipleCh0 << '\n';
    output << "ncal1_asymmetry_after_fraction=" << format_double(fraction(stats.ncal1Asymmetry)) << '\n';
    output << "lstilbene_asymmetry_after_fraction=" << format_double(fraction(stats.lStilbeneAsymmetry)) << '\n';
    output << "sstilbene_asymmetry_after_fraction=" << format_double(fraction(stats.sStilbeneAsymmetry)) << '\n';
    output << "elapsed_seconds=" << format_double(elapsedSeconds) << '\n';
}

void record_rate(std::unordered_map<std::int64_t, RateBinCounts>& bins,
                 AnalysisStats& stats,
                 std::int64_t absoluteHitNs,
                 int channel,
                 std::int64_t rateBinWidthNs) {
    const std::int64_t binIndex = floor_div(absoluteHitNs - stats.runStartTriggerNs, rateBinWidthNs);
    RateBinCounts& entry = bins[binIndex];
    if (channel == kCh0Channel) {
        entry.ch0 += 1;
    } else if (channel == kNcal1Channel) {
        entry.ncal1 += 1;
    } else if (channel == kLargeStilbeneChannel) {
        entry.lStilbene += 1;
    } else if (channel == kSmallStilbeneChannel) {
        entry.sStilbene += 1;
    }
    stats.minHitAbsoluteNs = std::min(stats.minHitAbsoluteNs, absoluteHitNs);
    stats.maxHitAbsoluteNs = std::max(stats.maxHitAbsoluteNs, absoluteHitNs);
    stats.minRateBinIndex = std::min(stats.minRateBinIndex, binIndex);
    stats.maxRateBinIndex = std::max(stats.maxRateBinIndex, binIndex);
}

void finalize_channel_window(const std::vector<double>& rises,
                             const std::vector<PulsePair>& pairs,
                             double referenceTimeNs,
                             const Config& config,
                             Histogram1D& foldedCounts,
                             Histogram2D& foldedTot,
                             WindowDetectorCounts& windowCounts) {
    for (const double riseNs : rises) {
        const double relativeNs = riseNs - referenceTimeNs;
        if (relativeNs < config.foldedMinNs || relativeNs >= config.foldedMaxNs) {
            continue;
        }
        fill_histogram(foldedCounts, relativeNs);
        windowCounts.total += 1;
        if (relativeNs >= 0.0) {
            windowCounts.after += 1;
        }
    }

    for (const PulsePair& pair : pairs) {
        const double relativeNs = pair.riseNs - referenceTimeNs;
        if (relativeNs < config.foldedMinNs || relativeNs >= config.foldedMaxNs) {
            continue;
        }
        if (pair.totNs < config.totMinNs || pair.totNs >= config.totMaxNs) {
            continue;
        }
        fill_histogram(foldedTot, relativeNs, pair.totNs);
    }
}

void finalize_signal_window(const WindowState& window,
                            std::uint64_t windowIndex,
                            const Config& config,
                            AnalysisStats& stats,
                            std::vector<AsymmetryBlock>& asymmetryBlocks,
                            Histogram1D& hNcal1Folded,
                            Histogram1D& hLStilbeneFolded,
                            Histogram1D& hSStilbeneFolded,
                            Histogram2D& hNcal1Tot,
                            Histogram2D& hLStilbeneTot,
                            Histogram2D& hSStilbeneTot) {
    stats.finalizedWindows += 1;

    AsymmetryBlock& block = ensure_block(asymmetryBlocks, windowIndex, config.asymmetryBlockSizeWindows);
    block.windowCount += 1;

    if (window.ch0Rises.empty()) {
        stats.windowsWithoutCh0 += 1;
        return;
    }

    stats.windowsWithValidCh0 += 1;
    block.validCh0Windows += 1;
    if (window.ch0Rises.size() > 1) {
        stats.windowsWithMultipleCh0 += 1;
    }

    const double referenceTimeNs = *std::min_element(window.ch0Rises.begin(), window.ch0Rises.end());
    WindowDetectorCounts ncal1Counts;
    WindowDetectorCounts lStilbeneCounts;
    WindowDetectorCounts sStilbeneCounts;

    finalize_channel_window(window.ncal1Rises, window.ncal1Pairs, referenceTimeNs, config, hNcal1Folded, hNcal1Tot, ncal1Counts);
    finalize_channel_window(window.lStilbeneRises, window.lStilbenePairs, referenceTimeNs, config, hLStilbeneFolded, hLStilbeneTot, lStilbeneCounts);
    finalize_channel_window(window.sStilbeneRises, window.sStilbenePairs, referenceTimeNs, config, hSStilbeneFolded, hSStilbeneTot, sStilbeneCounts);

    accumulate_counts(stats.ncal1Asymmetry, ncal1Counts);
    accumulate_counts(stats.lStilbeneAsymmetry, lStilbeneCounts);
    accumulate_counts(stats.sStilbeneAsymmetry, sStilbeneCounts);
    accumulate_counts(block.ncal1, ncal1Counts);
    accumulate_counts(block.lStilbene, lStilbeneCounts);
    accumulate_counts(block.sStilbene, sStilbeneCounts);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Config config = parse_args(argc, argv);
        std::ifstream input(config.inputPath);
        if (!input.is_open()) {
            throw std::runtime_error("Unable to open input file: " + config.inputPath);
        }

        Histogram1D hNcal1Folded = make_histogram_1d(
            "custom_RawRiseCountsNcal1FoldedCh0Ref",
            "Ncal1 raw rising-edge counts vs ch0-relative time",
            config.foldedBins,
            config.foldedMinNs,
            config.foldedMaxNs,
            "T relative to ch0 (ns)"
        );
        Histogram1D hLStilbeneFolded = make_histogram_1d(
            "custom_RawRiseCountsLstilbeneFoldedCh0Ref",
            "Lstilbene raw rising-edge counts vs ch0-relative time",
            config.foldedBins,
            config.foldedMinNs,
            config.foldedMaxNs,
            "T relative to ch0 (ns)"
        );
        Histogram1D hSStilbeneFolded = make_histogram_1d(
            "custom_RawRiseCountsSstilbeneFoldedCh0Ref",
            "Sstilbene raw rising-edge counts vs ch0-relative time",
            config.foldedBins,
            config.foldedMinNs,
            config.foldedMaxNs,
            "T relative to ch0 (ns)"
        );

        Histogram2D hNcal1Tot = make_histogram_2d(
            "custom_ProfNcal1Ch0Ref",
            "Ncal1 ToT vs ch0-relative time",
            config.foldedBins,
            config.foldedMinNs,
            config.foldedMaxNs,
            "T relative to ch0 (ns)",
            config.totBins,
            config.totMinNs,
            config.totMaxNs,
            "ToT (ns)"
        );
        Histogram2D hLStilbeneTot = make_histogram_2d(
            "custom_ProfLstilbeneCh0Ref",
            "Lstilbene ToT vs ch0-relative time",
            config.foldedBins,
            config.foldedMinNs,
            config.foldedMaxNs,
            "T relative to ch0 (ns)",
            config.totBins,
            config.totMinNs,
            config.totMaxNs,
            "ToT (ns)"
        );
        Histogram2D hSStilbeneTot = make_histogram_2d(
            "custom_ProfSstilbeneCh0Ref",
            "Sstilbene ToT vs ch0-relative time",
            config.foldedBins,
            config.foldedMinNs,
            config.foldedMaxNs,
            "T relative to ch0 (ns)",
            config.totBins,
            config.totMinNs,
            config.totMaxNs,
            "ToT (ns)"
        );

        AnalysisStats stats;
        std::unordered_map<std::int64_t, RateBinCounts> rateBins;
        std::vector<AsymmetryBlock> asymmetryBlocks;
        WindowState currentWindow;

        int currentTdcOrdinal = -1;
        bool haveLastTrigger0 = false;
        double lastTrigger0Seconds = 0.0;
        std::uint64_t currentWindowIndex = 0;
        bool haveSignalWindow = false;

        auto finalize_current_signal_window = [&]() {
            if (currentTdcOrdinal == kSignalOrdinal && haveSignalWindow) {
                finalize_signal_window(
                    currentWindow,
                    currentWindowIndex,
                    config,
                    stats,
                    asymmetryBlocks,
                    hNcal1Folded,
                    hLStilbeneFolded,
                    hSStilbeneFolded,
                    hNcal1Tot,
                    hLStilbeneTot,
                    hSStilbeneTot
                );
                currentWindow = WindowState{};
                currentWindowIndex += 1;
                haveSignalWindow = false;
            }
        };

        std::string line;
        const auto startedAt = std::chrono::steady_clock::now();
        while (std::getline(input, line)) {
            int tdcOrdinal = 0;
            int totalTdcs = 0;
            if (parse_header(line, tdcOrdinal, totalTdcs)) {
                finalize_current_signal_window();
                stats.parsedHeaders += 1;
                currentTdcOrdinal = tdcOrdinal;
                if (currentTdcOrdinal == kSignalOrdinal) {
                    haveSignalWindow = true;
                }
                continue;
            }

            if (const auto trigger = parse_trigger(line); trigger.has_value()) {
                if (currentTdcOrdinal == kTriggerOrdinal) {
                    if (haveLastTrigger0 && *trigger + config.triggerResetThresholdSeconds < lastTrigger0Seconds) {
                        stats.duplicateRestartDetected = true;
                        stats.duplicateRestartPreviousTriggerSeconds = lastTrigger0Seconds;
                        stats.duplicateRestartNewTriggerSeconds = *trigger;
                        stats.duplicateRestartHeadersSkipped = stats.parsedHeaders;
                        break;
                    }
                    haveLastTrigger0 = true;
                    lastTrigger0Seconds = *trigger;
                } else if (currentTdcOrdinal == kSignalOrdinal) {
                    currentWindow.hasTrigger = true;
                    currentWindow.triggerNs = static_cast<std::int64_t>(std::llround(*trigger * 1.0e9));
                    if (!stats.hasRunStartTrigger) {
                        stats.hasRunStartTrigger = true;
                        stats.runStartTriggerNs = currentWindow.triggerNs;
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
            if (currentTdcOrdinal != kSignalOrdinal || !haveSignalWindow) {
                continue;
            }
            if (channel != kCh0Channel && channel != kNcal1Channel && channel != kLargeStilbeneChannel && channel != kSmallStilbeneChannel) {
                continue;
            }
            if (!is_valid_hit_time_ns(timeNs)) {
                continue;
            }

            if (isRising == 1) {
                if (currentWindow.hasTrigger && stats.hasRunStartTrigger) {
                    const std::int64_t absoluteHitNs = currentWindow.triggerNs + static_cast<std::int64_t>(std::llround(timeNs));
                    record_rate(
                        rateBins,
                        stats,
                        absoluteHitNs,
                        channel,
                        static_cast<std::int64_t>(std::llround(config.rateBinWidthNs))
                    );
                }

                if (channel == kCh0Channel) {
                    currentWindow.ch0Rises.push_back(timeNs);
                    stats.validCh0Rises += 1;
                } else if (channel == kNcal1Channel) {
                    currentWindow.ncal1Rises.push_back(timeNs);
                    currentWindow.pendingNcal1Rises.push_back(timeNs);
                    stats.validNcal1Rises += 1;
                } else if (channel == kLargeStilbeneChannel) {
                    currentWindow.lStilbeneRises.push_back(timeNs);
                    currentWindow.pendingLStilbeneRises.push_back(timeNs);
                    stats.validLStilbeneRises += 1;
                } else if (channel == kSmallStilbeneChannel) {
                    currentWindow.sStilbeneRises.push_back(timeNs);
                    currentWindow.pendingSStilbeneRises.push_back(timeNs);
                    stats.validSStilbeneRises += 1;
                }
                continue;
            }

            auto consume_pair = [&](std::deque<double>& pending, std::vector<PulsePair>& pairs, std::uint64_t& validPairs) {
                if (pending.empty()) {
                    return;
                }
                const double riseNs = pending.front();
                pending.pop_front();
                const double totNs = timeNs - riseNs;
                if (!(totNs > 0.0) || !std::isfinite(totNs)) {
                    return;
                }
                pairs.push_back(PulsePair{riseNs, totNs});
                validPairs += 1;
            };

            if (channel == kNcal1Channel) {
                consume_pair(currentWindow.pendingNcal1Rises, currentWindow.ncal1Pairs, stats.validNcal1Pairs);
            } else if (channel == kLargeStilbeneChannel) {
                consume_pair(currentWindow.pendingLStilbeneRises, currentWindow.lStilbenePairs, stats.validLStilbenePairs);
            } else if (channel == kSmallStilbeneChannel) {
                consume_pair(currentWindow.pendingSStilbeneRises, currentWindow.sStilbenePairs, stats.validSStilbenePairs);
            }
        }

        finalize_current_signal_window();

        const auto endedAt = std::chrono::steady_clock::now();
        const double elapsedSeconds = std::chrono::duration_cast<std::chrono::duration<double>>(endedAt - startedAt).count();

        write_histogram_1d(config.outputPrefix + "_custom_RawRiseCountsNcal1FoldedCh0Ref_hist.txt", hNcal1Folded);
        write_histogram_1d(config.outputPrefix + "_custom_RawRiseCountsLstilbeneFoldedCh0Ref_hist.txt", hLStilbeneFolded);
        write_histogram_1d(config.outputPrefix + "_custom_RawRiseCountsSstilbeneFoldedCh0Ref_hist.txt", hSStilbeneFolded);
        write_histogram_2d(config.outputPrefix + "_custom_ProfNcal1Ch0Ref_hist.txt", hNcal1Tot);
        write_histogram_2d(config.outputPrefix + "_custom_ProfLstilbeneCh0Ref_hist.txt", hLStilbeneTot);
        write_histogram_2d(config.outputPrefix + "_custom_ProfSstilbeneCh0Ref_hist.txt", hSStilbeneTot);
        write_asymmetry_blocks_file(config.outputPrefix + "_asymmetry_blocks.txt", config, asymmetryBlocks);
        write_rates_file(config.outputPrefix + "_rates.txt", config, stats, rateBins);
        write_summary(config.outputPrefix + "_summary.txt", config, stats, asymmetryBlocks, elapsedSeconds);

        std::cout << "Finalized windows: " << stats.finalizedWindows << '\n';
        std::cout << "Asymmetry blocks: " << asymmetryBlocks.size() << '\n';
        std::cout << "Valid ch0 windows: " << stats.windowsWithValidCh0 << '\n';
        std::cout << "Ncal1 after fraction: "
                  << format_double(stats.ncal1Asymmetry.total == 0 ? 0.0 : static_cast<double>(stats.ncal1Asymmetry.after) / static_cast<double>(stats.ncal1Asymmetry.total))
                  << '\n';
        std::cout << "Lstilbene after fraction: "
                  << format_double(stats.lStilbeneAsymmetry.total == 0 ? 0.0 : static_cast<double>(stats.lStilbeneAsymmetry.after) / static_cast<double>(stats.lStilbeneAsymmetry.total))
                  << '\n';
        std::cout << "Sstilbene after fraction: "
                  << format_double(stats.sStilbeneAsymmetry.total == 0 ? 0.0 : static_cast<double>(stats.sStilbeneAsymmetry.after) / static_cast<double>(stats.sStilbeneAsymmetry.total))
                  << '\n';
        std::cout << "Duplicate restart detected: " << (stats.duplicateRestartDetected ? "true" : "false") << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}