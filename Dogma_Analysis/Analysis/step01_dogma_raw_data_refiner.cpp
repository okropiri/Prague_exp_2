#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace fs = std::filesystem;

namespace {

// Standalone preprocessing step for DOGMA raw text exports.
// It keeps only validated rise/fall pairs, writes a cleaned signal-window stream,
// and writes a pulse table that later analyses can read directly.

struct Config {
    std::string inputPath;
    std::string outputRoot;
    std::string runKey;
    std::string analysisSuffix = "rawRefined";
    std::optional<int> processedTdcOrdinal;
    std::optional<int> triggerTdcOrdinal;
    double triggerResetThresholdSeconds = 1.0;
    double maxAcceptedHitTimeNs = 50000.0;
    double totMinNs = 1.0;
    double totMaxNs = 5000.0;
    int maxAnomalySamplesPerGroup = 10;
    double progressReportIntervalSeconds = 15.0;
};

struct PendingRise {
    bool hasValue = false;
    int tdcIndexOneBased = 0;
    double riseNs = 0.0;
};

struct WindowState {
    std::string headerLine;
    int tdcOrdinal = -1;
    bool hasGlobalTriggerSeconds = false;
    double globalTriggerSeconds = 0.0;
    std::uint64_t windowIndex = 0;
    std::map<int, PendingRise> pendingRiseByChannel;
    std::vector<std::string> acceptedStreamRows;
};

struct ChannelStats {
    std::uint64_t totalRiseEntries = 0;
    std::uint64_t totalFallEntries = 0;
    std::uint64_t acceptedPairs = 0;
    std::uint64_t orphanedRises = 0;
    std::uint64_t orphanedFalls = 0;
    std::uint64_t replacedPendingRises = 0;
    std::uint64_t windowEndOrphanedRises = 0;
    std::uint64_t totBelowMin = 0;
    std::uint64_t totAboveMax = 0;
    std::uint64_t invalidHitTimeRows = 0;
    std::uint64_t unexpectedEdgeRows = 0;
};

struct RestartEvent {
    std::uint64_t occurrenceIndex = 0;
    std::uint64_t triggerLineNumber = 0;
    std::uint64_t parsedHeaderIndex = 0;
    std::uint64_t processedWindowIndex = 0;
    int tdcOrdinal = -1;
    double previousRawTriggerSeconds = 0.0;
    double newRawTriggerSeconds = 0.0;
    double previousAdjustedTriggerSeconds = 0.0;
    double newAdjustedTriggerSeconds = 0.0;
    double appliedOffsetSeconds = 0.0;
};

struct GlobalStats {
    std::uint64_t parsedHeaders = 0;
    std::uint64_t parsedRows = 0;
    std::uint64_t processedRows = 0;
    std::uint64_t processedWindowsSeen = 0;
    std::uint64_t finalizedWindows = 0;
    std::uint64_t refinedWindows = 0;
    bool duplicateRestartDetected = false;
    double duplicateRestartPreviousTriggerSeconds = 0.0;
    double duplicateRestartNewTriggerSeconds = 0.0;
    std::uint64_t duplicateRestartCount = 0;
    double duplicateRestartTimeOffsetSeconds = 0.0;
    std::uint64_t duplicateRestartHeadersSkipped = 0;
    std::vector<RestartEvent> restartEvents;
};

struct AnomalySample {
    std::uint64_t windowIndex = 0;
    int tdcOrdinal = -1;
    int channel = -1;
    std::string kind;
    bool hasGlobalTriggerSeconds = false;
    double globalTriggerSeconds = 0.0;
    int riseTdcIndexOneBased = 0;
    int fallTdcIndexOneBased = 0;
    double riseNs = std::numeric_limits<double>::quiet_NaN();
    double fallNs = std::numeric_limits<double>::quiet_NaN();
    double totNs = std::numeric_limits<double>::quiet_NaN();
    std::string detail;
};

bool starts_with(std::string_view text, std::string_view prefix) {
    return text.size() >= prefix.size() && text.substr(0, prefix.size()) == prefix;
}

bool ends_with(std::string_view text, std::string_view suffix) {
    return text.size() >= suffix.size() && text.substr(text.size() - suffix.size()) == suffix;
}

bool is_all_digits(std::string_view text) {
    if (text.empty()) {
        return false;
    }
    for (const char ch : text) {
        if (!std::isdigit(static_cast<unsigned char>(ch))) {
            return false;
        }
    }
    return true;
}

std::string format_double(double value) {
    if (!std::isfinite(value)) {
        return "nan";
    }
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(6) << value;
    std::string text = stream.str();
    while (!text.empty() && text.back() == '0') {
        text.pop_back();
    }
    if (!text.empty() && text.back() == '.') {
        text.pop_back();
    }
    if (text.empty()) {
        return "0";
    }
    return text;
}

std::string format_optional_int(const std::optional<int>& value) {
    return value.has_value() ? std::to_string(*value) : std::string("all");
}

std::string format_duration(double seconds) {
    if (!std::isfinite(seconds) || seconds < 0.0) {
        return "unknown";
    }
    const auto totalSeconds = static_cast<std::uint64_t>(std::llround(seconds));
    const std::uint64_t hours = totalSeconds / 3600;
    const std::uint64_t minutes = (totalSeconds % 3600) / 60;
    const std::uint64_t secs = totalSeconds % 60;
    std::ostringstream stream;
    stream << std::setfill('0') << std::setw(2) << hours << ':'
           << std::setw(2) << minutes << ':'
           << std::setw(2) << secs;
    return stream.str();
}

std::string format_bytes(std::uintmax_t bytes) {
    static constexpr const char* kUnits[] = {"B", "KB", "MB", "GB", "TB"};
    double value = static_cast<double>(bytes);
    std::size_t unitIndex = 0;
    while (value >= 1024.0 && unitIndex + 1 < std::size(kUnits)) {
        value /= 1024.0;
        unitIndex += 1;
    }
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(unitIndex == 0 ? 0 : 2) << value << ' ' << kUnits[unitIndex];
    return stream.str();
}

std::string format_fraction_percent(std::uint64_t numerator, std::uint64_t denominator) {
    if (denominator == 0) {
        return "0";
    }
    const double percent = 100.0 * static_cast<double>(numerator) / static_cast<double>(denominator);
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(3) << percent;
    return stream.str();
}

std::string strip_known_suffixes(std::string text) {
    bool changed = true;
    while (changed) {
        changed = false;
        if (ends_with(text, ".dld.dat")) {
            text.erase(text.size() - 8);
            changed = true;
        } else if (ends_with(text, ".dld")) {
            text.erase(text.size() - 4);
            changed = true;
        } else if (ends_with(text, ".dat")) {
            text.erase(text.size() - 4);
            changed = true;
        }
    }
    return text;
}

std::vector<std::string> split_tokens(const std::string& text) {
    std::vector<std::string> tokens;
    std::string current;
    for (const char ch : text) {
        if (ch == '_') {
            if (!current.empty()) {
                tokens.push_back(current);
                current.clear();
            }
            continue;
        }
        current.push_back(ch);
    }
    if (!current.empty()) {
        tokens.push_back(current);
    }
    return tokens;
}

std::string derive_run_key(const std::string& inputPath) {
    const fs::path path(inputPath);
    const std::string baseName = strip_known_suffixes(path.filename().string());
    const std::vector<std::string> tokens = split_tokens(baseName);
    if (tokens.empty()) {
        return "dogma_run";
    }

    const std::string detector = tokens.front();

    std::string windowToken;
    for (std::size_t index = 0; index < tokens.size(); ++index) {
        if (tokens[index] == "Pos") {
            break;
        }
        if (ends_with(tokens[index], "ns") || ends_with(tokens[index], "us") || ends_with(tokens[index], "ms")) {
            windowToken = tokens[index];
        }
    }

    std::string positionToken;
    for (std::size_t index = 0; index + 1 < tokens.size(); ++index) {
        if (tokens[index] == "Pos") {
            positionToken = "Pos_" + tokens[index + 1];
            break;
        }
    }

    std::string runId;
    for (std::size_t index = tokens.size(); index > 0; --index) {
        const std::string& token = tokens[index - 1];
        if (token.size() == 4 && is_all_digits(token)) {
            runId = token;
            break;
        }
    }

    if (!detector.empty() && !windowToken.empty() && !positionToken.empty() && !runId.empty()) {
        return detector + "_" + windowToken + "_" + positionToken + "_" + runId;
    }

    std::vector<std::string> fallbackTokens = tokens;
    if (fallbackTokens.size() >= 2 && fallbackTokens[fallbackTokens.size() - 2].size() == 8 &&
        is_all_digits(fallbackTokens[fallbackTokens.size() - 2]) &&
        fallbackTokens.back().size() == 4 && is_all_digits(fallbackTokens.back())) {
        fallbackTokens.erase(fallbackTokens.end() - 2);
    }

    std::ostringstream fallback;
    for (std::size_t index = 0; index < fallbackTokens.size(); ++index) {
        if (index != 0) {
            fallback << '_';
        }
        fallback << fallbackTokens[index];
    }
    return fallback.str();
}

void print_usage(const char* executableName) {
    std::cerr
        << "Usage:\n"
        << "  " << executableName << " --input <raw.dld.dat> --output-root <dir> [options]\n\n"
        << "What it writes:\n"
        << "  1. <run_key>_rawRefined/<run_key>_rawRefined.dld.dat\n"
        << "     Cleaned signal-window stream with only validated rise/fall pairs.\n"
        << "  2. <run_key>_rawRefined/<run_key>_rawRefined_pulses.tsv\n"
        << "     One accepted pulse per row for downstream analysis.\n"
        << "  3. <run_key>_rawRefined/<run_key>_rawRefined_validation.tsv\n"
        << "     Per-channel validation counts with an ALL row.\n"
        << "  4. <run_key>_rawRefined/<run_key>_rawRefined_anomaly_samples.tsv\n"
        << "     Small sample log of orphan and bad-ToT cases per channel.\n"
        << "  5. <run_key>_rawRefined/<run_key>_rawRefined_summary.txt\n"
        << "     Run-level metadata and output paths.\n\n"
        << "Main options:\n"
        << "  --run-key <text>                 Override the derived run key.\n"
        << "  --analysis-suffix <text>         Directory/file suffix. Default: rawRefined\n"
        << "  --signal-tdc-ordinal <int>       Only process this TDC ordinal. Default: all\n"
        << "  --trigger-tdc-ordinal <int>      Only use trigger lines from this TDC ordinal. Default: all\n"
        << "  --tot-min-ns <value>             Minimum accepted ToT, inclusive. Default: 1\n"
        << "  --tot-max-ns <value>             Maximum accepted ToT, inclusive. Default: 5000\n"
        << "  --max-anomaly-samples <int>      Keep this many samples per channel and anomaly kind. Default: 10\n"
        << "  --progress-report-interval-seconds <v>  Progress print interval. Default: 15\n"
        << "  --max-accepted-hit-time-ns <v>   Reject non-finite or larger |time| rows. Default: 50000\n"
        << "  --trigger-reset-threshold-seconds <v>  Restart rollback threshold. Default: 1\n";
}

Config parse_args(int argc, char** argv) {
    Config config;
    for (int index = 1; index < argc; ++index) {
        const std::string arg = argv[index];
        if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        }
        if (arg == "--input" && index + 1 < argc) {
            config.inputPath = argv[++index];
        } else if (arg == "--output-root" && index + 1 < argc) {
            config.outputRoot = argv[++index];
        } else if (arg == "--run-key" && index + 1 < argc) {
            config.runKey = argv[++index];
        } else if (arg == "--analysis-suffix" && index + 1 < argc) {
            config.analysisSuffix = argv[++index];
        } else if (arg == "--signal-tdc-ordinal" && index + 1 < argc) {
            config.processedTdcOrdinal = std::stoi(argv[++index]);
        } else if (arg == "--trigger-tdc-ordinal" && index + 1 < argc) {
            config.triggerTdcOrdinal = std::stoi(argv[++index]);
        } else if (arg == "--tot-min-ns" && index + 1 < argc) {
            config.totMinNs = std::stod(argv[++index]);
        } else if (arg == "--tot-max-ns" && index + 1 < argc) {
            config.totMaxNs = std::stod(argv[++index]);
        } else if (arg == "--max-anomaly-samples" && index + 1 < argc) {
            config.maxAnomalySamplesPerGroup = std::stoi(argv[++index]);
        } else if (arg == "--progress-report-interval-seconds" && index + 1 < argc) {
            config.progressReportIntervalSeconds = std::stod(argv[++index]);
        } else if (arg == "--max-accepted-hit-time-ns" && index + 1 < argc) {
            config.maxAcceptedHitTimeNs = std::stod(argv[++index]);
        } else if (arg == "--trigger-reset-threshold-seconds" && index + 1 < argc) {
            config.triggerResetThresholdSeconds = std::stod(argv[++index]);
        } else {
            throw std::runtime_error("Unknown or incomplete argument: " + arg);
        }
    }

    if (config.inputPath.empty()) {
        throw std::runtime_error("Missing required argument: --input");
    }
    if (config.outputRoot.empty()) {
        throw std::runtime_error("Missing required argument: --output-root");
    }
    if (config.totMinNs < 0.0) {
        throw std::runtime_error("--tot-min-ns must be non-negative");
    }
    if (config.totMaxNs < config.totMinNs) {
        throw std::runtime_error("--tot-max-ns must be >= --tot-min-ns");
    }
    if (config.maxAnomalySamplesPerGroup < 0) {
        throw std::runtime_error("--max-anomaly-samples must be >= 0");
    }
    if (config.progressReportIntervalSeconds <= 0.0) {
        throw std::runtime_error("--progress-report-interval-seconds must be > 0");
    }
    if (config.runKey.empty()) {
        config.runKey = derive_run_key(config.inputPath);
    }
    return config;
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

bool parse_row(const std::string& line, int& tdcIndexOneBased, int& channel, int& edge, double& timeNs) {
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
    edge = static_cast<int>(third);
    timeNs = fourth;
    return true;
}

bool is_valid_hit_time_ns(double timeNs, double maxAcceptedHitTimeNs) {
    return std::isfinite(timeNs) && std::abs(timeNs) <= maxAcceptedHitTimeNs;
}

bool should_track_tdc(const std::optional<int>& filter, int tdcOrdinal) {
    return !filter.has_value() || tdcOrdinal == *filter;
}

std::string make_anomaly_group_key(int channel, const std::string& kind) {
    return std::to_string(channel) + "|" + kind;
}

void record_anomaly_sample(std::vector<AnomalySample>& anomalySamples,
                           std::map<std::string, int>& anomalySampleCounts,
                           const Config& config,
                           const WindowState& window,
                           int channel,
                           const std::string& kind,
                           int riseTdcIndexOneBased,
                           int fallTdcIndexOneBased,
                           double riseNs,
                           double fallNs,
                           double totNs,
                           const std::string& detail) {
    const std::string key = make_anomaly_group_key(channel, kind);
    int& count = anomalySampleCounts[key];
    if (count >= config.maxAnomalySamplesPerGroup) {
        return;
    }
    anomalySamples.push_back(AnomalySample{
        window.windowIndex,
        window.tdcOrdinal,
        channel,
        kind,
        window.hasGlobalTriggerSeconds,
        window.globalTriggerSeconds,
        riseTdcIndexOneBased,
        fallTdcIndexOneBased,
        riseNs,
        fallNs,
        totNs,
        detail,
    });
    count += 1;
}

ChannelStats sum_channel_stats(const std::map<int, ChannelStats>& statsByChannel) {
    ChannelStats total;
    for (const auto& [channel, stats] : statsByChannel) {
        (void)channel;
        total.totalRiseEntries += stats.totalRiseEntries;
        total.totalFallEntries += stats.totalFallEntries;
        total.acceptedPairs += stats.acceptedPairs;
        total.orphanedRises += stats.orphanedRises;
        total.orphanedFalls += stats.orphanedFalls;
        total.replacedPendingRises += stats.replacedPendingRises;
        total.windowEndOrphanedRises += stats.windowEndOrphanedRises;
        total.totBelowMin += stats.totBelowMin;
        total.totAboveMax += stats.totAboveMax;
        total.invalidHitTimeRows += stats.invalidHitTimeRows;
        total.unexpectedEdgeRows += stats.unexpectedEdgeRows;
    }
    return total;
}

void write_pulse_table_header(std::ofstream& output, const Config& config) {
    output << "# input_file=" << config.inputPath << '\n';
    output << "# run_key=" << config.runKey << '\n';
    output << "# tot_min_ns=" << format_double(config.totMinNs) << '\n';
    output << "# tot_max_ns=" << format_double(config.totMaxNs) << '\n';
    output << "window_index\ttdc_ordinal\tglobal_trigger_seconds\tchannel\trise_tdc_index_one_based\tfall_tdc_index_one_based\trise_ns\tfall_ns\ttot_ns\n";
}

void write_pulse_record(std::ofstream& output,
                        std::uint64_t windowIndex,
                        int tdcOrdinal,
                        bool hasGlobalTriggerSeconds,
                        double globalTriggerSeconds,
                        int channel,
                        int riseTdcIndexOneBased,
                        int fallTdcIndexOneBased,
                        double riseNs,
                        double fallNs,
                        double totNs) {
    output << windowIndex << '\t'
           << tdcOrdinal << '\t';
    if (hasGlobalTriggerSeconds) {
        output << format_double(globalTriggerSeconds);
    }
    output << '\t'
           << channel << '\t'
           << riseTdcIndexOneBased << '\t'
           << fallTdcIndexOneBased << '\t'
           << format_double(riseNs) << '\t'
           << format_double(fallNs) << '\t'
           << format_double(totNs) << '\n';
}

void write_validation_table(const fs::path& outputPath,
                            const Config& config,
                            const GlobalStats& globalStats,
                            const std::map<int, ChannelStats>& statsByChannel) {
    std::ofstream output(outputPath);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write validation table: " + outputPath.string());
    }

    const ChannelStats total = sum_channel_stats(statsByChannel);

    output << "# input_file=" << config.inputPath << '\n';
    output << "# run_key=" << config.runKey << '\n';
    output << "# processed_tdc_ordinal=" << format_optional_int(config.processedTdcOrdinal) << '\n';
    output << "# trigger_tdc_ordinal=" << format_optional_int(config.triggerTdcOrdinal) << '\n';
    output << "# tot_min_ns=" << format_double(config.totMinNs) << '\n';
    output << "# tot_max_ns=" << format_double(config.totMaxNs) << '\n';
    output << "# max_accepted_hit_time_ns=" << format_double(config.maxAcceptedHitTimeNs) << '\n';
    output << "# duplicate_restart_detected=" << (globalStats.duplicateRestartDetected ? "true" : "false") << '\n';
    output << "channel\ttot_min_ns\ttot_max_ns\ttotal_rise_entries\ttotal_fall_entries\taccepted_pairs\torphaned_rises\torphaned_rises_percent_of_total_rises\torphaned_falls\torphaned_falls_percent_of_total_falls\ttot_below_min\ttot_above_max\tmisvalidated_tots\tinvalid_hit_time_rows\treplaced_pending_rises\twindow_end_orphaned_rises\tunexpected_edge_rows\n";

    const auto write_row = [&](const std::string& channelLabel, const ChannelStats& stats) {
        output << channelLabel << '\t'
               << format_double(config.totMinNs) << '\t'
               << format_double(config.totMaxNs) << '\t'
               << stats.totalRiseEntries << '\t'
               << stats.totalFallEntries << '\t'
               << stats.acceptedPairs << '\t'
               << stats.orphanedRises << '\t'
               << format_fraction_percent(stats.orphanedRises, stats.totalRiseEntries) << '\t'
               << stats.orphanedFalls << '\t'
               << format_fraction_percent(stats.orphanedFalls, stats.totalFallEntries) << '\t'
               << stats.totBelowMin << '\t'
               << stats.totAboveMax << '\t'
               << (stats.totBelowMin + stats.totAboveMax) << '\t'
               << stats.invalidHitTimeRows << '\t'
               << stats.replacedPendingRises << '\t'
               << stats.windowEndOrphanedRises << '\t'
               << stats.unexpectedEdgeRows << '\n';
    };

    write_row("ALL", total);
    for (const auto& [channel, stats] : statsByChannel) {
        write_row(std::to_string(channel), stats);
    }
}

void write_anomaly_samples(const fs::path& outputPath,
                           const Config& config,
                           const std::vector<AnomalySample>& anomalySamples) {
    std::ofstream output(outputPath);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write anomaly sample table: " + outputPath.string());
    }

    output << "# input_file=" << config.inputPath << '\n';
    output << "# run_key=" << config.runKey << '\n';
    output << "# processed_tdc_ordinal=" << format_optional_int(config.processedTdcOrdinal) << '\n';
    output << "# trigger_tdc_ordinal=" << format_optional_int(config.triggerTdcOrdinal) << '\n';
    output << "# tot_min_ns=" << format_double(config.totMinNs) << '\n';
    output << "# tot_max_ns=" << format_double(config.totMaxNs) << '\n';
    output << "# max_anomaly_samples_per_group=" << config.maxAnomalySamplesPerGroup << '\n';
    output << "window_index\ttdc_ordinal\tchannel\tanomaly_kind\tglobal_trigger_seconds\trise_tdc_index_one_based\tfall_tdc_index_one_based\trise_ns\tfall_ns\ttot_ns\tdetail\n";

    for (const AnomalySample& sample : anomalySamples) {
        output << sample.windowIndex << '\t'
               << sample.tdcOrdinal << '\t'
               << sample.channel << '\t'
               << sample.kind << '\t';
        if (sample.hasGlobalTriggerSeconds) {
            output << format_double(sample.globalTriggerSeconds);
        }
        output << '\t'
               << sample.riseTdcIndexOneBased << '\t'
               << sample.fallTdcIndexOneBased << '\t'
               << format_double(sample.riseNs) << '\t'
               << format_double(sample.fallNs) << '\t'
               << format_double(sample.totNs) << '\t'
               << sample.detail << '\n';
    }
}

void write_summary(const fs::path& outputPath,
                   const Config& config,
                   const GlobalStats& globalStats,
                   const std::map<int, ChannelStats>& statsByChannel,
                   const fs::path& refinedStreamPath,
                   const fs::path& pulseTablePath,
                   const fs::path& validationTablePath,
                   const fs::path& anomalyTablePath,
                   std::uintmax_t inputBytesTotal,
                   double elapsedSeconds) {
    std::ofstream output(outputPath);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write summary file: " + outputPath.string());
    }

    const ChannelStats total = sum_channel_stats(statsByChannel);

    output << "input_file=" << config.inputPath << '\n';
    output << "run_key=" << config.runKey << '\n';
    output << "processed_tdc_ordinal=" << format_optional_int(config.processedTdcOrdinal) << '\n';
    output << "trigger_tdc_ordinal=" << format_optional_int(config.triggerTdcOrdinal) << '\n';
    output << "tot_min_ns=" << format_double(config.totMinNs) << '\n';
    output << "tot_max_ns=" << format_double(config.totMaxNs) << '\n';
    output << "max_accepted_hit_time_ns=" << format_double(config.maxAcceptedHitTimeNs) << '\n';
    output << "progress_report_interval_seconds=" << format_double(config.progressReportIntervalSeconds) << '\n';
    output << "input_size_bytes=" << inputBytesTotal << '\n';
    output << "parsed_headers=" << globalStats.parsedHeaders << '\n';
    output << "parsed_rows=" << globalStats.parsedRows << '\n';
    output << "processed_rows=" << globalStats.processedRows << '\n';
    output << "processed_windows_seen=" << globalStats.processedWindowsSeen << '\n';
    output << "finalized_windows=" << globalStats.finalizedWindows << '\n';
    output << "refined_windows=" << globalStats.refinedWindows << '\n';
    output << "duplicate_restart_detected=" << (globalStats.duplicateRestartDetected ? "true" : "false") << '\n';
    output << "duplicate_restart_previous_trigger_seconds=" << format_double(globalStats.duplicateRestartPreviousTriggerSeconds) << '\n';
    output << "duplicate_restart_new_trigger_seconds=" << format_double(globalStats.duplicateRestartNewTriggerSeconds) << '\n';
    output << "duplicate_restart_count=" << globalStats.duplicateRestartCount << '\n';
    output << "duplicate_restart_time_offset_seconds=" << format_double(globalStats.duplicateRestartTimeOffsetSeconds) << '\n';
    output << "duplicate_restart_headers_skipped=" << globalStats.duplicateRestartHeadersSkipped << '\n';
    for (const RestartEvent& event : globalStats.restartEvents) {
        const std::string prefix = "restart_event_" + std::to_string(event.occurrenceIndex);
        output << prefix << "_trigger_line_number=" << event.triggerLineNumber << '\n';
        output << prefix << "_parsed_header_index=" << event.parsedHeaderIndex << '\n';
        output << prefix << "_processed_window_index=" << event.processedWindowIndex << '\n';
        output << prefix << "_tdc_ordinal=" << event.tdcOrdinal << '\n';
        output << prefix << "_previous_raw_trigger_seconds=" << format_double(event.previousRawTriggerSeconds) << '\n';
        output << prefix << "_new_raw_trigger_seconds=" << format_double(event.newRawTriggerSeconds) << '\n';
        output << prefix << "_previous_adjusted_trigger_seconds=" << format_double(event.previousAdjustedTriggerSeconds) << '\n';
        output << prefix << "_new_adjusted_trigger_seconds=" << format_double(event.newAdjustedTriggerSeconds) << '\n';
        output << prefix << "_applied_offset_seconds=" << format_double(event.appliedOffsetSeconds) << '\n';
    }
    output << "total_rise_entries=" << total.totalRiseEntries << '\n';
    output << "total_fall_entries=" << total.totalFallEntries << '\n';
    output << "accepted_pairs=" << total.acceptedPairs << '\n';
    output << "orphaned_rises=" << total.orphanedRises << '\n';
    output << "orphaned_falls=" << total.orphanedFalls << '\n';
    output << "tot_below_min=" << total.totBelowMin << '\n';
    output << "tot_above_max=" << total.totAboveMax << '\n';
    output << "misvalidated_tots=" << (total.totBelowMin + total.totAboveMax) << '\n';
    output << "invalid_hit_time_rows=" << total.invalidHitTimeRows << '\n';
    output << "orphaned_rises_percent_of_total_rises=" << format_fraction_percent(total.orphanedRises, total.totalRiseEntries) << '\n';
    output << "orphaned_falls_percent_of_total_falls=" << format_fraction_percent(total.orphanedFalls, total.totalFallEntries) << '\n';
    output << "elapsed_seconds=" << format_double(elapsedSeconds) << '\n';
    output << "refined_stream_output=" << refinedStreamPath.string() << '\n';
    output << "pulse_table_output=" << pulseTablePath.string() << '\n';
    output << "validation_table_output=" << validationTablePath.string() << '\n';
    output << "anomaly_samples_output=" << anomalyTablePath.string() << '\n';

    output << "\n# How to read this summary\n";
    output << "# This file starts with machine-readable key=value lines.\n";
    output << "# The notes below explain the most important fields in plain language.\n";
    output << "#\n";
    output << "# Example of one accepted pulse\n";
    output << "#   rise row: 2 21 1 272.267\n";
    output << "#   fall row: 2 21 0 281.278\n";
    output << "#   ToT = 281.278 - 272.267 = 9.011 ns\n";
    output << "#   This contributes one accepted pair on channel 21.\n";
    output << "#\n";
    output << "# Field guide\n";
    output << "# input_file: source DOGMA text stream that was cleaned.\n";
    output << "# processed_tdc_ordinal: which TDC ordinal was cleaned. 'all' means every TDC block was processed.\n";
    output << "# trigger_tdc_ordinal: which TDC ordinal supplied GlobalTriggerTime values. 'all' means any block can supply them.\n";
    output << "# tot_min_ns / tot_max_ns: accepted ToT window. Pairs outside this window are rejected.\n";
    output << "# max_accepted_hit_time_ns: sanity limit for a single hit time. Larger absolute values are treated as bad rows.\n";
    output << "# parsed_headers: count of DOGMA block headers starting with 'TDC ...'.\n";
    output << "# parsed_rows: count of raw four-column data rows that were parsed successfully.\n";
    output << "# processed_rows: parsed rows that belonged to the selected TDC blocks.\n";
    output << "# processed_windows_seen: number of processed DOGMA blocks.\n";
    output << "# finalized_windows: processed blocks that were closed and checked by the cleaner.\n";
    output << "# refined_windows: processed blocks that produced at least one accepted pair in the cleaned output stream.\n";
    output << "# duplicate_restart_detected: true means GlobalTriggerTime jumped backwards by more than the restart threshold.\n";
    output << "# duplicate_restart_count: how many backward trigger jumps were detected and stitched forward.\n";
    output << "# duplicate_restart_time_offset_seconds: cumulative offset added to later raw trigger values so cleaned output keeps a monotonic absolute-time axis.\n";
    output << "# duplicate_restart_headers_skipped: always 0 in stitched mode because cleaning continues after a restart.\n";
    output << "# restart_event_<n>_*: exact details for each detected restart, including the trigger line number in the raw text file and the raw trigger values before and after the jump.\n";
    output << "# total_rise_entries: all rise rows seen in processed blocks, before ToT validation.\n";
    output << "# total_fall_entries: all fall rows seen in processed blocks, before ToT validation.\n";
    output << "# accepted_pairs: rise-fall pairs that passed the sanity checks and ToT window.\n";
    output << "# orphaned_rises: rise rows that never became an accepted pair.\n";
    output << "# orphaned_falls: fall rows that arrived while no rise was open on that channel.\n";
    output << "# tot_below_min: pairs rejected because the computed ToT was below the configured minimum.\n";
    output << "# tot_above_max: pairs rejected because the computed ToT was above the configured maximum.\n";
    output << "# misvalidated_tots: total ToT rejections = tot_below_min + tot_above_max.\n";
    output << "# invalid_hit_time_rows: rows rejected before pairing because the time value was not sane.\n";
    output << "# orphaned_rises_percent_of_total_rises: orphaned_rises divided by total_rise_entries, in percent.\n";
    output << "# orphaned_falls_percent_of_total_falls: orphaned_falls divided by total_fall_entries, in percent.\n";
    output << "# elapsed_seconds: total wall-clock time spent in this cleaning run.\n";
    output << "#\n";
    output << "# Breakdown of orphaned rises\n";
    output << "# replaced_pending_rises: a newer rise arrived before any fall, so the older rise was discarded.\n";
    output << "# window_end_orphaned_rises: a rise was still open when the DOGMA block ended.\n";
    output << "# In the current logic: orphaned_rises = replaced_pending_rises + window_end_orphaned_rises.\n";
    if (!globalStats.restartEvents.empty()) {
        output << "#\n";
        output << "# Restart events\n";
        output << "# occurrence trigger_line parsed_header processed_window tdc previous_raw_s new_raw_s previous_adjusted_s new_adjusted_s applied_offset_s\n";
        for (const RestartEvent& event : globalStats.restartEvents) {
            output << "# " << event.occurrenceIndex << ' '
                   << event.triggerLineNumber << ' '
                   << event.parsedHeaderIndex << ' '
                   << event.processedWindowIndex << ' '
                   << event.tdcOrdinal << ' '
                   << format_double(event.previousRawTriggerSeconds) << ' '
                   << format_double(event.newRawTriggerSeconds) << ' '
                   << format_double(event.previousAdjustedTriggerSeconds) << ' '
                   << format_double(event.newAdjustedTriggerSeconds) << ' '
                   << format_double(event.appliedOffsetSeconds) << '\n';
        }
    }
    output << "#\n";
    output << "# Output files\n";
    output << "# refined_stream_output: cleaned DOGMA-style text stream.\n";
    output << "# pulse_table_output: one accepted pulse per row, easier for downstream analysis.\n";
    output << "# validation_table_output: per-channel statistics table with one ALL row, including orphan percentages next to the orphan columns.\n";
    output << "# anomaly_samples_output: a small sample of representative bad rows and rejected pairs.\n";
}

void finalize_window(WindowState& window,
                     GlobalStats& globalStats,
                     std::map<int, ChannelStats>& statsByChannel,
                     std::vector<AnomalySample>& anomalySamples,
                     std::map<std::string, int>& anomalySampleCounts,
                     const Config& config,
                     std::ofstream& refinedStreamOutput) {
    globalStats.finalizedWindows += 1;

    for (auto& [channel, pendingRise] : window.pendingRiseByChannel) {
        if (!pendingRise.hasValue) {
            continue;
        }
        ChannelStats& stats = statsByChannel[channel];
        stats.orphanedRises += 1;
        stats.windowEndOrphanedRises += 1;
        record_anomaly_sample(anomalySamples,
                              anomalySampleCounts,
                              config,
                              window,
                              channel,
                              "orphan_rise_window_end",
                              pendingRise.tdcIndexOneBased,
                              0,
                              pendingRise.riseNs,
                              std::numeric_limits<double>::quiet_NaN(),
                              std::numeric_limits<double>::quiet_NaN(),
                              "pending rise remained open at end of block");
        pendingRise = PendingRise{};
    }

    if (!window.acceptedStreamRows.empty()) {
        refinedStreamOutput << window.headerLine << '\n';
        if (window.hasGlobalTriggerSeconds) {
            refinedStreamOutput << "GlobalTriggerTime " << format_double(window.globalTriggerSeconds) << '\n';
        }
        for (const std::string& row : window.acceptedStreamRows) {
            refinedStreamOutput << row << '\n';
        }
        globalStats.refinedWindows += 1;
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Config config = parse_args(argc, argv);
        const auto startedAt = std::chrono::steady_clock::now();

        std::ifstream input(config.inputPath);
        if (!input.is_open()) {
            throw std::runtime_error("Unable to open input file: " + config.inputPath);
        }

        std::error_code fileSizeError;
        const std::uintmax_t inputBytesTotal = fs::file_size(config.inputPath, fileSizeError);
        const bool haveInputSize = !fileSizeError;

        const fs::path outputDir = fs::path(config.outputRoot) / (config.runKey + "_" + config.analysisSuffix);
        fs::create_directories(outputDir);
        const std::string prefix = config.runKey + "_" + config.analysisSuffix;
        const fs::path refinedStreamPath = outputDir / (prefix + ".dld.dat");
        const fs::path pulseTablePath = outputDir / (prefix + "_pulses.tsv");
        const fs::path validationTablePath = outputDir / (prefix + "_validation.tsv");
        const fs::path anomalyTablePath = outputDir / (prefix + "_anomaly_samples.tsv");
        const fs::path summaryPath = outputDir / (prefix + "_summary.txt");

        std::ofstream refinedStreamOutput(refinedStreamPath);
        if (!refinedStreamOutput.is_open()) {
            throw std::runtime_error("Unable to open refined stream output: " + refinedStreamPath.string());
        }

        std::ofstream pulseTableOutput(pulseTablePath);
        if (!pulseTableOutput.is_open()) {
            throw std::runtime_error("Unable to open pulse table output: " + pulseTablePath.string());
        }
        write_pulse_table_header(pulseTableOutput, config);

        GlobalStats globalStats;
        std::map<int, ChannelStats> statsByChannel;
        std::vector<AnomalySample> anomalySamples;
        std::map<std::string, int> anomalySampleCounts;
        WindowState currentWindow;

        int currentTdcOrdinal = -1;
        bool haveLastTrigger = false;
        double lastRawTriggerSeconds = 0.0;
        double lastAdjustedTriggerSeconds = 0.0;
        double triggerTimeOffsetSeconds = 0.0;
        bool haveProcessedWindow = false;
        auto lastProgressAt = startedAt;
        std::uint64_t progressCheckCounter = 0;

        auto finalize_current_window = [&]() {
            if (haveProcessedWindow) {
                finalize_window(currentWindow,
                                globalStats,
                                statsByChannel,
                                anomalySamples,
                                anomalySampleCounts,
                                config,
                                refinedStreamOutput);
                currentWindow = WindowState{};
                haveProcessedWindow = false;
            }
        };

        auto report_progress = [&](bool force) {
            const auto now = std::chrono::steady_clock::now();
            const double elapsedSeconds = std::chrono::duration_cast<std::chrono::duration<double>>(now - startedAt).count();
            const double sinceLastReport = std::chrono::duration_cast<std::chrono::duration<double>>(now - lastProgressAt).count();
            if (!force && sinceLastReport < config.progressReportIntervalSeconds) {
                return;
            }
            lastProgressAt = now;

            std::uintmax_t bytesRead = 0;
            const std::streampos tellPosition = input.tellg();
            if (tellPosition != std::streampos(-1)) {
                bytesRead = static_cast<std::uintmax_t>(tellPosition);
            } else if (haveInputSize) {
                bytesRead = inputBytesTotal;
            }

            const ChannelStats total = sum_channel_stats(statsByChannel);
            std::cerr << "[progress] ";
            if (haveInputSize && inputBytesTotal > 0) {
                const double percent = 100.0 * static_cast<double>(bytesRead) / static_cast<double>(inputBytesTotal);
                const double bytesPerSecond = elapsedSeconds > 0.0 ? static_cast<double>(bytesRead) / elapsedSeconds : 0.0;
                const double etaSeconds = (bytesPerSecond > 0.0 && bytesRead < inputBytesTotal)
                    ? static_cast<double>(inputBytesTotal - bytesRead) / bytesPerSecond
                    : 0.0;
                std::cerr << std::fixed << std::setprecision(2)
                          << percent << "% | "
                          << format_bytes(bytesRead) << "/" << format_bytes(inputBytesTotal)
                          << std::defaultfloat;
                if (etaSeconds > 0.0) {
                    std::cerr << " | eta " << format_duration(etaSeconds);
                }
            } else {
                std::cerr << format_bytes(bytesRead) << " read";
            }
            std::cerr << " | elapsed " << format_duration(elapsedSeconds)
                      << " | windows " << globalStats.processedWindowsSeen
                      << " | pairs " << total.acceptedPairs
                      << " | orphan_rises " << total.orphanedRises
                      << " | orphan_falls " << total.orphanedFalls
                      << '\n';
        };

        std::cout << "Input size: " << (haveInputSize ? format_bytes(inputBytesTotal) : std::string("unknown")) << '\n';
        std::cout << "Progress interval: " << format_double(config.progressReportIntervalSeconds) << " s" << '\n';

        std::string line;
        std::uint64_t currentLineNumber = 0;
        while (std::getline(input, line)) {
            currentLineNumber += 1;
            int tdcOrdinal = 0;
            int totalTdcs = 0;
            if (parse_header(line, tdcOrdinal, totalTdcs)) {
                finalize_current_window();
                (void)totalTdcs;
                globalStats.parsedHeaders += 1;
                currentTdcOrdinal = tdcOrdinal;
                if (should_track_tdc(config.processedTdcOrdinal, currentTdcOrdinal)) {
                    haveProcessedWindow = true;
                    currentWindow = WindowState{};
                    currentWindow.headerLine = line;
                    currentWindow.tdcOrdinal = currentTdcOrdinal;
                    currentWindow.windowIndex = globalStats.processedWindowsSeen + 1;
                    globalStats.processedWindowsSeen += 1;
                    if (haveLastTrigger) {
                        currentWindow.hasGlobalTriggerSeconds = true;
                        currentWindow.globalTriggerSeconds = lastAdjustedTriggerSeconds;
                    }
                }
                continue;
            }

            if (const auto trigger = parse_trigger(line); trigger.has_value()) {
                if (should_track_tdc(config.triggerTdcOrdinal, currentTdcOrdinal)) {
                    double adjustedTriggerSeconds = *trigger + triggerTimeOffsetSeconds;
                    if (haveLastTrigger && *trigger + config.triggerResetThresholdSeconds < lastRawTriggerSeconds) {
                        globalStats.duplicateRestartDetected = true;
                        globalStats.duplicateRestartPreviousTriggerSeconds = lastRawTriggerSeconds;
                        globalStats.duplicateRestartNewTriggerSeconds = *trigger;
                        globalStats.duplicateRestartCount += 1;
                        const double previousAdjustedTriggerSeconds = lastAdjustedTriggerSeconds;
                        triggerTimeOffsetSeconds = lastAdjustedTriggerSeconds - *trigger;
                        adjustedTriggerSeconds = *trigger + triggerTimeOffsetSeconds;
                        globalStats.duplicateRestartTimeOffsetSeconds = triggerTimeOffsetSeconds;
                        globalStats.duplicateRestartHeadersSkipped = 0;
                        globalStats.restartEvents.push_back(RestartEvent{
                            globalStats.duplicateRestartCount,
                            currentLineNumber,
                            globalStats.parsedHeaders,
                            globalStats.processedWindowsSeen,
                            currentTdcOrdinal,
                            lastRawTriggerSeconds,
                            *trigger,
                            previousAdjustedTriggerSeconds,
                            adjustedTriggerSeconds,
                            triggerTimeOffsetSeconds,
                        });
                    }
                    haveLastTrigger = true;
                    lastRawTriggerSeconds = *trigger;
                    lastAdjustedTriggerSeconds = adjustedTriggerSeconds;
                    if (haveProcessedWindow) {
                        currentWindow.hasGlobalTriggerSeconds = true;
                        currentWindow.globalTriggerSeconds = adjustedTriggerSeconds;
                    }
                }
                continue;
            }

            int tdcIndexOneBased = 0;
            int channel = 0;
            int edge = 0;
            double timeNs = 0.0;
            if (!parse_row(line, tdcIndexOneBased, channel, edge, timeNs)) {
                continue;
            }

            globalStats.parsedRows += 1;
            progressCheckCounter += 1;
            if ((progressCheckCounter % 250000) == 0) {
                report_progress(false);
            }
            if (!haveProcessedWindow) {
                continue;
            }

            globalStats.processedRows += 1;
            ChannelStats& stats = statsByChannel[channel];

            if (edge == 1) {
                stats.totalRiseEntries += 1;
            } else if (edge == 0) {
                stats.totalFallEntries += 1;
            } else {
                stats.unexpectedEdgeRows += 1;
                continue;
            }

            if (!is_valid_hit_time_ns(timeNs, config.maxAcceptedHitTimeNs)) {
                stats.invalidHitTimeRows += 1;
                record_anomaly_sample(anomalySamples,
                                      anomalySampleCounts,
                                      config,
                                      currentWindow,
                                      channel,
                                      "invalid_hit_time",
                                      edge == 1 ? tdcIndexOneBased : 0,
                                      edge == 0 ? tdcIndexOneBased : 0,
                                      edge == 1 ? timeNs : std::numeric_limits<double>::quiet_NaN(),
                                      edge == 0 ? timeNs : std::numeric_limits<double>::quiet_NaN(),
                                      std::numeric_limits<double>::quiet_NaN(),
                                      "row rejected by hit-time sanity limit");
                continue;
            }

            PendingRise& pendingRise = currentWindow.pendingRiseByChannel[channel];
            // Keep one open pulse per channel. If a second rise appears before a fall,
            // the older rise is treated as orphaned and the new rise becomes the only candidate.
            if (edge == 1) {
                if (pendingRise.hasValue) {
                    stats.orphanedRises += 1;
                    stats.replacedPendingRises += 1;
                    record_anomaly_sample(anomalySamples,
                                          anomalySampleCounts,
                                          config,
                                          currentWindow,
                                          channel,
                                          "orphan_rise_replaced",
                                          pendingRise.tdcIndexOneBased,
                                          tdcIndexOneBased,
                                          pendingRise.riseNs,
                                          timeNs,
                                          std::numeric_limits<double>::quiet_NaN(),
                                          "new rise replaced an older open rise");
                }
                pendingRise = PendingRise{true, tdcIndexOneBased, timeNs};
                continue;
            }

            if (!pendingRise.hasValue) {
                stats.orphanedFalls += 1;
                record_anomaly_sample(anomalySamples,
                                      anomalySampleCounts,
                                      config,
                                      currentWindow,
                                      channel,
                                      "orphan_fall",
                                      0,
                                      tdcIndexOneBased,
                                      std::numeric_limits<double>::quiet_NaN(),
                                      timeNs,
                                      std::numeric_limits<double>::quiet_NaN(),
                                      "fall arrived with no open rise");
                continue;
            }

            const PendingRise rise = pendingRise;
            pendingRise = PendingRise{};
            const double totNs = timeNs - rise.riseNs;
            if (!std::isfinite(totNs) || totNs < config.totMinNs) {
                stats.totBelowMin += 1;
                record_anomaly_sample(anomalySamples,
                                      anomalySampleCounts,
                                      config,
                                      currentWindow,
                                      channel,
                                      "tot_below_min",
                                      rise.tdcIndexOneBased,
                                      tdcIndexOneBased,
                                      rise.riseNs,
                                      timeNs,
                                      totNs,
                                      "pair rejected because ToT is below the configured minimum");
                continue;
            }
            if (totNs > config.totMaxNs) {
                stats.totAboveMax += 1;
                record_anomaly_sample(anomalySamples,
                                      anomalySampleCounts,
                                      config,
                                      currentWindow,
                                      channel,
                                      "tot_above_max",
                                      rise.tdcIndexOneBased,
                                      tdcIndexOneBased,
                                      rise.riseNs,
                                      timeNs,
                                      totNs,
                                      "pair rejected because ToT is above the configured maximum");
                continue;
            }

            stats.acceptedPairs += 1;
            write_pulse_record(pulseTableOutput,
                               currentWindow.windowIndex,
                               currentWindow.tdcOrdinal,
                               currentWindow.hasGlobalTriggerSeconds,
                               currentWindow.globalTriggerSeconds,
                               channel,
                               rise.tdcIndexOneBased,
                               tdcIndexOneBased,
                               rise.riseNs,
                               timeNs,
                               totNs);

            currentWindow.acceptedStreamRows.push_back(
                std::to_string(rise.tdcIndexOneBased) + " " + std::to_string(channel) + " 1 " + format_double(rise.riseNs));
            currentWindow.acceptedStreamRows.push_back(
                std::to_string(tdcIndexOneBased) + " " + std::to_string(channel) + " 0 " + format_double(timeNs));
        }

        finalize_current_window();
        report_progress(true);

        const auto endedAt = std::chrono::steady_clock::now();
        const double elapsedSeconds = std::chrono::duration_cast<std::chrono::duration<double>>(endedAt - startedAt).count();

        write_validation_table(validationTablePath, config, globalStats, statsByChannel);
        write_anomaly_samples(anomalyTablePath, config, anomalySamples);
        write_summary(summaryPath,
                      config,
                      globalStats,
                      statsByChannel,
                      refinedStreamPath,
                      pulseTablePath,
                      validationTablePath,
                      anomalyTablePath,
                      haveInputSize ? inputBytesTotal : 0,
                      elapsedSeconds);

        std::cout << "Refined stream: " << refinedStreamPath << '\n';
        std::cout << "Pulse table: " << pulseTablePath << '\n';
        std::cout << "Validation table: " << validationTablePath << '\n';
        std::cout << "Anomaly samples: " << anomalyTablePath << '\n';
        std::cout << "Summary: " << summaryPath << '\n';
        std::cout << "Run key: " << config.runKey << '\n';
        std::cout << "Elapsed: " << format_duration(elapsedSeconds) << '\n';
        std::cout << "Trigger restarts detected: " << globalStats.duplicateRestartCount << '\n';
        for (const RestartEvent& event : globalStats.restartEvents) {
            std::cout << "  restart " << event.occurrenceIndex
                      << " at line " << event.triggerLineNumber
                      << ": raw " << format_double(event.previousRawTriggerSeconds)
                      << " -> " << format_double(event.newRawTriggerSeconds)
                      << " s, adjusted " << format_double(event.previousAdjustedTriggerSeconds)
                      << " -> " << format_double(event.newAdjustedTriggerSeconds)
                      << " s, offset " << format_double(event.appliedOffsetSeconds)
                      << " s" << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}