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
    std::uint64_t blockSizeWindows = 1000;
    double triggerResetThresholdSeconds = 1.0;
    double windowMinNs = -5000.0;
    double windowMaxNs = 5000.0;
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

struct WindowCounts {
    WindowDetectorCounts ncal1;
    WindowDetectorCounts lStilbene;
    WindowDetectorCounts sStilbene;
};

struct BlockAggregate {
    std::uint64_t windowCount = 0;
    AggregateCounts ncal1;
    AggregateCounts lStilbene;
    AggregateCounts sStilbene;
    AggregateCounts lStilbeneMult1;
    AggregateCounts lStilbeneMult2To3;
    AggregateCounts lStilbeneMult4To7;
    AggregateCounts lStilbeneMult8Plus;
    std::uint64_t lStilbeneZeroWindows = 0;
    std::uint64_t lStilbeneMult1Windows = 0;
    std::uint64_t lStilbeneMult2To3Windows = 0;
    std::uint64_t lStilbeneMult4To7Windows = 0;
    std::uint64_t lStilbeneMult8PlusWindows = 0;
};

struct AnalysisStats {
    std::uint64_t parsedHeaders = 0;
    std::uint64_t parsedRows = 0;
    std::uint64_t finalizedWindows = 0;
    bool duplicateRestartDetected = false;
    double duplicateRestartPreviousTriggerSeconds = 0.0;
    double duplicateRestartNewTriggerSeconds = 0.0;
    std::uint64_t duplicateRestartHeadersSkipped = 0;
    AggregateCounts ncal1;
    AggregateCounts lStilbene;
    AggregateCounts sStilbene;
    AggregateCounts lStilbeneMult1;
    AggregateCounts lStilbeneMult2To3;
    AggregateCounts lStilbeneMult4To7;
    AggregateCounts lStilbeneMult8Plus;
    std::uint64_t lStilbeneZeroWindows = 0;
    std::uint64_t lStilbeneMult1Windows = 0;
    std::uint64_t lStilbeneMult2To3Windows = 0;
    std::uint64_t lStilbeneMult4To7Windows = 0;
    std::uint64_t lStilbeneMult8PlusWindows = 0;
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
        } else if (arg == "--block-size-windows" && index + 1 < argc) {
            config.blockSizeWindows = static_cast<std::uint64_t>(std::stoull(argv[++index]));
        } else if (arg == "--trigger-reset-threshold-seconds" && index + 1 < argc) {
            config.triggerResetThresholdSeconds = std::stod(argv[++index]);
        } else if (arg == "--window-min-ns" && index + 1 < argc) {
            config.windowMinNs = std::stod(argv[++index]);
        } else if (arg == "--window-max-ns" && index + 1 < argc) {
            config.windowMaxNs = std::stod(argv[++index]);
        } else if (arg == "--help") {
            std::cout
                << "Usage: dogma_trgref_asymmetry_investigation --input <file> --output-prefix <prefix>\n"
                << "  [--block-size-windows 1000] [--trigger-reset-threshold-seconds 1.0]\n"
                << "  [--window-min-ns -5000] [--window-max-ns 5000]\n";
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
    if (config.blockSizeWindows == 0) {
        throw std::runtime_error("block-size-windows must be positive");
    }
    if (config.triggerResetThresholdSeconds <= 0.0) {
        throw std::runtime_error("trigger-reset-threshold-seconds must be positive");
    }
    if (!(config.windowMinNs < config.windowMaxNs)) {
        throw std::runtime_error("window-min-ns must be smaller than window-max-ns");
    }
    return config;
}

std::string format_double(double value) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(6) << value;
    return stream.str();
}

double safe_ratio(std::uint64_t numerator, std::uint64_t denominator) {
    if (denominator == 0) {
        return 0.0;
    }
    return static_cast<double>(numerator) / static_cast<double>(denominator);
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

bool is_in_asymmetry_window(double timeNs, const Config& config) {
    return timeNs >= config.windowMinNs && timeNs < config.windowMaxNs;
}

void accumulate_counts(AggregateCounts& aggregate, const WindowDetectorCounts& counts) {
    aggregate.total += counts.total;
    aggregate.after += counts.after;
    if (counts.total > 0) {
        aggregate.nonzeroWindows += 1;
    }
}

BlockAggregate& ensure_block(std::vector<BlockAggregate>& blocks,
                             std::uint64_t windowIndex,
                             std::uint64_t blockSizeWindows) {
    const std::size_t blockIndex = static_cast<std::size_t>(windowIndex / blockSizeWindows);
    if (blocks.size() <= blockIndex) {
        blocks.resize(blockIndex + 1);
    }
    return blocks[blockIndex];
}

void finalize_window(const WindowCounts& windowCounts,
                     std::uint64_t windowIndex,
                     const Config& config,
                     AnalysisStats& stats,
                     std::vector<BlockAggregate>& blocks) {
    stats.finalizedWindows += 1;

    BlockAggregate& block = ensure_block(blocks, windowIndex, config.blockSizeWindows);
    block.windowCount += 1;

    accumulate_counts(stats.ncal1, windowCounts.ncal1);
    accumulate_counts(stats.lStilbene, windowCounts.lStilbene);
    accumulate_counts(stats.sStilbene, windowCounts.sStilbene);

    accumulate_counts(block.ncal1, windowCounts.ncal1);
    accumulate_counts(block.lStilbene, windowCounts.lStilbene);
    accumulate_counts(block.sStilbene, windowCounts.sStilbene);

    const std::uint32_t lMultiplicity = windowCounts.lStilbene.total;
    if (lMultiplicity == 0) {
        stats.lStilbeneZeroWindows += 1;
        block.lStilbeneZeroWindows += 1;
    } else if (lMultiplicity == 1) {
        stats.lStilbeneMult1Windows += 1;
        block.lStilbeneMult1Windows += 1;
        accumulate_counts(stats.lStilbeneMult1, windowCounts.lStilbene);
        accumulate_counts(block.lStilbeneMult1, windowCounts.lStilbene);
    } else if (lMultiplicity <= 3) {
        stats.lStilbeneMult2To3Windows += 1;
        block.lStilbeneMult2To3Windows += 1;
        accumulate_counts(stats.lStilbeneMult2To3, windowCounts.lStilbene);
        accumulate_counts(block.lStilbeneMult2To3, windowCounts.lStilbene);
    } else if (lMultiplicity <= 7) {
        stats.lStilbeneMult4To7Windows += 1;
        block.lStilbeneMult4To7Windows += 1;
        accumulate_counts(stats.lStilbeneMult4To7, windowCounts.lStilbene);
        accumulate_counts(block.lStilbeneMult4To7, windowCounts.lStilbene);
    } else {
        stats.lStilbeneMult8PlusWindows += 1;
        block.lStilbeneMult8PlusWindows += 1;
        accumulate_counts(stats.lStilbeneMult8Plus, windowCounts.lStilbene);
        accumulate_counts(block.lStilbeneMult8Plus, windowCounts.lStilbene);
    }
}

void record_row(WindowDetectorCounts& counts, double timeNs) {
    counts.total += 1;
    if (timeNs >= 0.0) {
        counts.after += 1;
    }
}

void write_blocks_file(const std::string& path,
                       const Config& config,
                       const std::vector<BlockAggregate>& blocks) {
    std::ofstream output(path);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write block file: " + path);
    }
    output << "# block_size_windows=" << config.blockSizeWindows << '\n';
    output << "# window_min_ns=" << format_double(config.windowMinNs) << '\n';
    output << "# window_max_ns=" << format_double(config.windowMaxNs) << '\n';
    output << "# columns: block_index window_start window_end_exclusive windows_in_block"
           << " ncal1_after ncal1_total lstilbene_after lstilbene_total sstilbene_after sstilbene_total"
           << " lstilbene_mult1_after lstilbene_mult1_total"
           << " lstilbene_mult2to3_after lstilbene_mult2to3_total"
           << " lstilbene_mult4to7_after lstilbene_mult4to7_total"
           << " lstilbene_mult8plus_after lstilbene_mult8plus_total"
           << " lstilbene_zero_windows lstilbene_mult1_windows lstilbene_mult2to3_windows"
           << " lstilbene_mult4to7_windows lstilbene_mult8plus_windows\n";

    for (std::size_t blockIndex = 0; blockIndex < blocks.size(); ++blockIndex) {
        const BlockAggregate& block = blocks[blockIndex];
        const std::uint64_t windowStart = static_cast<std::uint64_t>(blockIndex) * config.blockSizeWindows;
        const std::uint64_t windowEndExclusive = windowStart + block.windowCount;
        output << blockIndex << ' '
               << windowStart << ' '
               << windowEndExclusive << ' '
               << block.windowCount << ' '
               << block.ncal1.after << ' ' << block.ncal1.total << ' '
               << block.lStilbene.after << ' ' << block.lStilbene.total << ' '
               << block.sStilbene.after << ' ' << block.sStilbene.total << ' '
               << block.lStilbeneMult1.after << ' ' << block.lStilbeneMult1.total << ' '
               << block.lStilbeneMult2To3.after << ' ' << block.lStilbeneMult2To3.total << ' '
               << block.lStilbeneMult4To7.after << ' ' << block.lStilbeneMult4To7.total << ' '
               << block.lStilbeneMult8Plus.after << ' ' << block.lStilbeneMult8Plus.total << ' '
               << block.lStilbeneZeroWindows << ' '
               << block.lStilbeneMult1Windows << ' '
               << block.lStilbeneMult2To3Windows << ' '
               << block.lStilbeneMult4To7Windows << ' '
               << block.lStilbeneMult8PlusWindows << '\n';
    }
}

void write_summary(const std::string& path,
                   const Config& config,
                   const AnalysisStats& stats,
                   std::size_t blockCount,
                   double elapsedSeconds) {
    std::ofstream output(path);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write summary file: " + path);
    }
    output << "input_file=" << config.inputPath << '\n';
    output << "block_size_windows=" << config.blockSizeWindows << '\n';
    output << "window_min_ns=" << format_double(config.windowMinNs) << '\n';
    output << "window_max_ns=" << format_double(config.windowMaxNs) << '\n';
    output << "trigger_reset_threshold_seconds=" << format_double(config.triggerResetThresholdSeconds) << '\n';
    output << "parsed_headers=" << stats.parsedHeaders << '\n';
    output << "parsed_rows=" << stats.parsedRows << '\n';
    output << "finalized_windows=" << stats.finalizedWindows << '\n';
    output << "block_count=" << blockCount << '\n';
    output << "duplicate_restart_detected=" << (stats.duplicateRestartDetected ? "true" : "false") << '\n';
    output << "duplicate_restart_previous_trigger_seconds=" << format_double(stats.duplicateRestartPreviousTriggerSeconds) << '\n';
    output << "duplicate_restart_new_trigger_seconds=" << format_double(stats.duplicateRestartNewTriggerSeconds) << '\n';
    output << "duplicate_restart_headers_skipped=" << stats.duplicateRestartHeadersSkipped << '\n';

    output << "ncal1_total_counts=" << stats.ncal1.total << '\n';
    output << "ncal1_after_counts=" << stats.ncal1.after << '\n';
    output << "ncal1_after_fraction=" << format_double(safe_ratio(stats.ncal1.after, stats.ncal1.total)) << '\n';
    output << "ncal1_nonzero_windows=" << stats.ncal1.nonzeroWindows << '\n';

    output << "lstilbene_total_counts=" << stats.lStilbene.total << '\n';
    output << "lstilbene_after_counts=" << stats.lStilbene.after << '\n';
    output << "lstilbene_after_fraction=" << format_double(safe_ratio(stats.lStilbene.after, stats.lStilbene.total)) << '\n';
    output << "lstilbene_nonzero_windows=" << stats.lStilbene.nonzeroWindows << '\n';

    output << "sstilbene_total_counts=" << stats.sStilbene.total << '\n';
    output << "sstilbene_after_counts=" << stats.sStilbene.after << '\n';
    output << "sstilbene_after_fraction=" << format_double(safe_ratio(stats.sStilbene.after, stats.sStilbene.total)) << '\n';
    output << "sstilbene_nonzero_windows=" << stats.sStilbene.nonzeroWindows << '\n';

    output << "lstilbene_zero_windows=" << stats.lStilbeneZeroWindows << '\n';
    output << "lstilbene_mult1_windows=" << stats.lStilbeneMult1Windows << '\n';
    output << "lstilbene_mult2to3_windows=" << stats.lStilbeneMult2To3Windows << '\n';
    output << "lstilbene_mult4to7_windows=" << stats.lStilbeneMult4To7Windows << '\n';
    output << "lstilbene_mult8plus_windows=" << stats.lStilbeneMult8PlusWindows << '\n';

    output << "lstilbene_mult1_total_counts=" << stats.lStilbeneMult1.total << '\n';
    output << "lstilbene_mult1_after_counts=" << stats.lStilbeneMult1.after << '\n';
    output << "lstilbene_mult1_after_fraction=" << format_double(safe_ratio(stats.lStilbeneMult1.after, stats.lStilbeneMult1.total)) << '\n';
    output << "lstilbene_mult2to3_total_counts=" << stats.lStilbeneMult2To3.total << '\n';
    output << "lstilbene_mult2to3_after_counts=" << stats.lStilbeneMult2To3.after << '\n';
    output << "lstilbene_mult2to3_after_fraction=" << format_double(safe_ratio(stats.lStilbeneMult2To3.after, stats.lStilbeneMult2To3.total)) << '\n';
    output << "lstilbene_mult4to7_total_counts=" << stats.lStilbeneMult4To7.total << '\n';
    output << "lstilbene_mult4to7_after_counts=" << stats.lStilbeneMult4To7.after << '\n';
    output << "lstilbene_mult4to7_after_fraction=" << format_double(safe_ratio(stats.lStilbeneMult4To7.after, stats.lStilbeneMult4To7.total)) << '\n';
    output << "lstilbene_mult8plus_total_counts=" << stats.lStilbeneMult8Plus.total << '\n';
    output << "lstilbene_mult8plus_after_counts=" << stats.lStilbeneMult8Plus.after << '\n';
    output << "lstilbene_mult8plus_after_fraction=" << format_double(safe_ratio(stats.lStilbeneMult8Plus.after, stats.lStilbeneMult8Plus.total)) << '\n';
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
        std::vector<BlockAggregate> blocks;
        WindowCounts currentWindowCounts;
        bool windowActive = false;
        bool stopOnDuplicate = false;
        int currentTdcOrdinal = -1;
        std::uint64_t currentWindowIndex = 0;
        bool haveLastWindowTrigger = false;
        double lastWindowTriggerSeconds = 0.0;
        std::string line;
        const auto startedAt = std::chrono::steady_clock::now();

        while (std::getline(input, line)) {
            int tdcOrdinal = 0;
            int totalTdcs = 0;
            if (parse_header(line, tdcOrdinal, totalTdcs)) {
                stats.parsedHeaders += 1;
                if (tdcOrdinal == 0) {
                    if (windowActive) {
                        finalize_window(currentWindowCounts, currentWindowIndex, config, stats, blocks);
                        currentWindowCounts = WindowCounts{};
                        currentWindowIndex += 1;
                    } else {
                        windowActive = true;
                        currentWindowIndex = 0;
                    }
                }
                currentTdcOrdinal = tdcOrdinal;
                continue;
            }

            if (const auto trigger = parse_trigger(line); trigger.has_value()) {
                if (currentTdcOrdinal == 0) {
                    if (haveLastWindowTrigger &&
                        *trigger + config.triggerResetThresholdSeconds < lastWindowTriggerSeconds) {
                        stats.duplicateRestartDetected = true;
                        stats.duplicateRestartPreviousTriggerSeconds = lastWindowTriggerSeconds;
                        stats.duplicateRestartNewTriggerSeconds = *trigger;
                        stats.duplicateRestartHeadersSkipped = stats.parsedHeaders;
                        stopOnDuplicate = true;
                        break;
                    }
                    haveLastWindowTrigger = true;
                    lastWindowTriggerSeconds = *trigger;
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
            if (!windowActive || currentTdcOrdinal != kTargetTdcOrdinal) {
                continue;
            }
            if (isRising != 1) {
                continue;
            }
            if (!is_valid_hit_time_ns(timeNs) || !is_in_asymmetry_window(timeNs, config)) {
                continue;
            }

            if (channel == kNcal1Channel) {
                record_row(currentWindowCounts.ncal1, timeNs);
            } else if (channel == kLargeStilbeneChannel) {
                record_row(currentWindowCounts.lStilbene, timeNs);
            } else if (channel == kSmallStilbeneChannel) {
                record_row(currentWindowCounts.sStilbene, timeNs);
            }
        }

        if (windowActive && !stopOnDuplicate) {
            finalize_window(currentWindowCounts, currentWindowIndex, config, stats, blocks);
        }

        const auto endedAt = std::chrono::steady_clock::now();
        const double elapsedSeconds = std::chrono::duration_cast<std::chrono::duration<double>>(endedAt - startedAt).count();

        write_blocks_file(config.outputPrefix + "_blocks.txt", config, blocks);
        write_summary(config.outputPrefix + "_summary.txt", config, stats, blocks.size(), elapsedSeconds);

        std::cout << "Finalized windows: " << stats.finalizedWindows << '\n';
        std::cout << "Block count: " << blocks.size() << '\n';
        std::cout << "Ncal1 after fraction: " << format_double(safe_ratio(stats.ncal1.after, stats.ncal1.total)) << '\n';
        std::cout << "Lstilbene after fraction: " << format_double(safe_ratio(stats.lStilbene.after, stats.lStilbene.total)) << '\n';
        std::cout << "Sstilbene after fraction: " << format_double(safe_ratio(stats.sStilbene.after, stats.sStilbene.total)) << '\n';
        std::cout << "Duplicate restart detected: " << (stats.duplicateRestartDetected ? "true" : "false") << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}