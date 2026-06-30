#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace {

constexpr double kSentinelValue = 2748779069440.0;
constexpr std::size_t kSampleLimit = 12;
constexpr std::size_t kResolvedHistoryLimit = 256;

struct Config {
    std::string inputPath;
    std::string summaryOutputPath;
    std::string histogramOutputPath;
    std::string deltaProfileOutputPath;
    std::string periodSummaryOutputPath;
    int referenceChannel = 0;
    int signalChannel = 2;
    double deltaMinNs = 0.0;
    double deltaMaxNs = 45000.0;
    double maxParentDeltaNs = -1.0;
    int deltaBins = 1125;
    double totMinNs = 0.0;
    double totMaxNs = 128.0;
    int totBins = 512;
    double profileBinWidthNs = 0.5;
    double profileTotMinNs = 0.0;
    double profileTotMaxNs = 128.0;
    double detrendWindowNs = 400.0;
    double periodSearchMinNs = 35.0;
    double periodSearchMaxNs = 45.0;
};

struct PeriodCandidate {
    double periodNs = 0.0;
    double score = 0.0;
    int lagBins = 0;
};

struct PeriodEstimate {
    bool available = false;
    double bestPeriodNs = 0.0;
    double bestScore = 0.0;
    std::uint64_t profiledPulses = 0;
    std::vector<double> detrendedCounts;
    std::vector<PeriodCandidate> topCandidates;
};

struct PulseStats {
    std::uint64_t rows = 0;
    std::uint64_t rises = 0;
    std::uint64_t falls = 0;
    std::uint64_t sentinelRows = 0;
    std::uint64_t duplicateRise = 0;
    std::uint64_t orphanRiseEndBlock = 0;
    std::uint64_t orphanFallWithoutRise = 0;
    std::uint64_t negativeTot = 0;
    std::uint64_t zeroTot = 0;
    std::uint64_t validPositivePairs = 0;
    double validTotMin = std::numeric_limits<double>::infinity();
    double validTotMax = -std::numeric_limits<double>::infinity();
    double validTotSum = 0.0;
};

struct AcceptedStats {
    std::uint64_t count = 0;
    double deltaMin = std::numeric_limits<double>::infinity();
    double deltaMax = -std::numeric_limits<double>::infinity();
    double deltaSum = 0.0;
    double totMin = std::numeric_limits<double>::infinity();
    double totMax = -std::numeric_limits<double>::infinity();
    double totSum = 0.0;
};

struct PendingPulse {
    bool hasValue = false;
    double triggerSeconds = 0.0;
    double leadingNs = 0.0;
    double absoluteLeadingNs = 0.0;
    int firstCol = 0;
};

struct ResolvedParent {
    bool exists = false;
    bool usable = false;
    double triggerSeconds = 0.0;
    double leadingNs = 0.0;
    double absoluteLeadingNs = 0.0;
    double totNs = 0.0;
    int firstCol = 0;
    std::string reason;
};

struct SignalPulse {
    double triggerSeconds = 0.0;
    double leadingNs = 0.0;
    double absoluteLeadingNs = 0.0;
    double totNs = 0.0;
    int firstCol = 0;
};

struct StreamState {
    PendingPulse pendingReference;
    PendingPulse pendingSignal;
    std::vector<ResolvedParent> resolvedReferenceHistory;
    std::vector<SignalPulse> waitingSignalPulses;
};

struct AnalysisStats {
    PulseStats reference;
    PulseStats signal;
    AcceptedStats accepted;
    std::uint64_t completedBlocks = 0;
    std::uint64_t acceptedSignalPulses = 0;
    std::uint64_t rejectedNoEarlierParent = 0;
    std::uint64_t rejectedBadParent = 0;
    std::uint64_t rejectedParentTooOld = 0;
    std::uint64_t rejectedNegativeParentDelta = 0;
    std::uint64_t histogrammedSignalPulses = 0;
    std::uint64_t deltaOverflow = 0;
    std::uint64_t histogramTotOverflow = 0;
    std::map<std::string, std::uint64_t> acceptedByStream;
    std::map<std::string, std::uint64_t> rejectedByReason;
    std::vector<std::string> sampleReferenceAnomalies;
    std::vector<std::string> sampleSignalAnomalies;
    std::vector<std::string> sampleParentRejects;
};

bool starts_with(std::string_view text, std::string_view prefix) {
    return text.size() >= prefix.size() && text.substr(0, prefix.size()) == prefix;
}

Config parse_args(int argc, char** argv) {
    Config config;
    for (int index = 1; index < argc; ++index) {
        const std::string arg = argv[index];
        if (arg == "--input" && index + 1 < argc) {
            config.inputPath = argv[++index];
        } else if (arg == "--summary-output" && index + 1 < argc) {
            config.summaryOutputPath = argv[++index];
        } else if (arg == "--histogram-output" && index + 1 < argc) {
            config.histogramOutputPath = argv[++index];
        } else if (arg == "--delta-profile-output" && index + 1 < argc) {
            config.deltaProfileOutputPath = argv[++index];
        } else if (arg == "--period-summary-output" && index + 1 < argc) {
            config.periodSummaryOutputPath = argv[++index];
        } else if (arg == "--reference-channel" && index + 1 < argc) {
            config.referenceChannel = std::stoi(argv[++index]);
        } else if (arg == "--signal-channel" && index + 1 < argc) {
            config.signalChannel = std::stoi(argv[++index]);
        } else if (arg == "--delta-min-ns" && index + 1 < argc) {
            config.deltaMinNs = std::stod(argv[++index]);
        } else if (arg == "--delta-max-ns" && index + 1 < argc) {
            config.deltaMaxNs = std::stod(argv[++index]);
        } else if (arg == "--max-parent-delta-ns" && index + 1 < argc) {
            config.maxParentDeltaNs = std::stod(argv[++index]);
        } else if (arg == "--delta-bins" && index + 1 < argc) {
            config.deltaBins = std::stoi(argv[++index]);
        } else if (arg == "--tot-min-ns" && index + 1 < argc) {
            config.totMinNs = std::stod(argv[++index]);
        } else if (arg == "--tot-max-ns" && index + 1 < argc) {
            config.totMaxNs = std::stod(argv[++index]);
        } else if (arg == "--tot-bins" && index + 1 < argc) {
            config.totBins = std::stoi(argv[++index]);
        } else if (arg == "--profile-bin-width-ns" && index + 1 < argc) {
            config.profileBinWidthNs = std::stod(argv[++index]);
        } else if (arg == "--profile-tot-min-ns" && index + 1 < argc) {
            config.profileTotMinNs = std::stod(argv[++index]);
        } else if (arg == "--profile-tot-max-ns" && index + 1 < argc) {
            config.profileTotMaxNs = std::stod(argv[++index]);
        } else if (arg == "--detrend-window-ns" && index + 1 < argc) {
            config.detrendWindowNs = std::stod(argv[++index]);
        } else if (arg == "--period-search-min-ns" && index + 1 < argc) {
            config.periodSearchMinNs = std::stod(argv[++index]);
        } else if (arg == "--period-search-max-ns" && index + 1 < argc) {
            config.periodSearchMaxNs = std::stod(argv[++index]);
        } else if (arg == "--help") {
            std::cout
                << "Usage: dogma_ch2_parent_time_tot_no_modulo --input <file> --summary-output <txt> --histogram-output <dat>\n"
                << "  [--reference-channel 0] [--signal-channel 2]\n"
                << "  [--delta-min-ns 0] [--delta-max-ns 45000] [--max-parent-delta-ns delta-max-ns] [--delta-bins 1125]\n"
                << "  [--tot-min-ns 0] [--tot-max-ns 128] [--tot-bins 512]\n"
                << "  [--delta-profile-output <tsv>] [--period-summary-output <txt>]\n"
                << "  [--profile-bin-width-ns 0.5] [--profile-tot-min-ns 0] [--profile-tot-max-ns 128]\n"
                << "  [--detrend-window-ns 400] [--period-search-min-ns 35] [--period-search-max-ns 45]\n";
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown or incomplete argument: " + arg);
        }
    }

    if (config.inputPath.empty()) {
        throw std::runtime_error("Missing required argument: --input");
    }
    if (config.summaryOutputPath.empty()) {
        throw std::runtime_error("Missing required argument: --summary-output");
    }
    if (config.histogramOutputPath.empty()) {
        throw std::runtime_error("Missing required argument: --histogram-output");
    }
    if (config.deltaBins <= 0 || config.totBins <= 0) {
        throw std::runtime_error("Histogram bin counts must be positive");
    }
    if (config.deltaMaxNs <= config.deltaMinNs) {
        throw std::runtime_error("delta-max-ns must be larger than delta-min-ns");
    }
    if (config.maxParentDeltaNs <= 0.0) {
        config.maxParentDeltaNs = config.deltaMaxNs;
    }
    if (config.maxParentDeltaNs <= 0.0) {
        throw std::runtime_error("max-parent-delta-ns must be positive");
    }
    if (config.totMaxNs <= config.totMinNs) {
        throw std::runtime_error("tot-max-ns must be larger than tot-min-ns");
    }
    if (config.profileBinWidthNs <= 0.0) {
        throw std::runtime_error("profile-bin-width-ns must be positive");
    }
    if (config.profileTotMaxNs <= config.profileTotMinNs) {
        throw std::runtime_error("profile-tot-max-ns must be larger than profile-tot-min-ns");
    }
    if (config.detrendWindowNs <= 0.0) {
        throw std::runtime_error("detrend-window-ns must be positive");
    }
    if (config.periodSearchMaxNs <= config.periodSearchMinNs) {
        throw std::runtime_error("period-search-max-ns must be larger than period-search-min-ns");
    }
    return config;
}

std::string parse_tdc_id(const std::string& line) {
    const std::string key = "TDC TDC_";
    const auto start = line.find(key);
    if (start == std::string::npos) {
        return "UNKNOWN";
    }
    const auto idStart = start + 4;
    auto idEnd = line.find(' ', idStart);
    if (idEnd == std::string::npos) {
        idEnd = line.size();
    }
    return line.substr(idStart, idEnd - idStart);
}

std::optional<double> parse_trigger(const std::string& line) {
    constexpr std::string_view prefix = "GlobalTriggerTime ";
    if (!starts_with(line, prefix)) {
        return std::nullopt;
    }
    const char* first = line.c_str() + prefix.size();
    char* end = nullptr;
    const double value = std::strtod(first, &end);
    return end == first ? std::nullopt : std::optional<double>(value);
}

bool parse_row(const std::string& line, int& firstCol, int& channel, int& edge, double& value) {
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

    firstCol = static_cast<int>(first);
    channel = static_cast<int>(second);
    edge = static_cast<int>(third);
    value = fourth;
    return true;
}

double absolute_time_ns(double triggerSeconds, double localNs) {
    return triggerSeconds * 1.0e9 + localNs;
}

std::string format_double(double value) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(6) << value;
    return stream.str();
}

void add_sample(std::vector<std::string>& samples, const std::string& text) {
    if (samples.size() < kSampleLimit) {
        samples.push_back(text);
    }
}

void update_positive_tot(PulseStats& stats, double totNs) {
    stats.validPositivePairs += 1;
    stats.validTotMin = std::min(stats.validTotMin, totNs);
    stats.validTotMax = std::max(stats.validTotMax, totNs);
    stats.validTotSum += totNs;
}

PeriodEstimate estimate_period(const Config& config,
                               const std::vector<std::uint64_t>& profileCounts) {
    PeriodEstimate result;
    if (profileCounts.empty()) {
        return result;
    }

    result.detrendedCounts.resize(profileCounts.size(), 0.0);
    std::vector<double> prefix(profileCounts.size() + 1, 0.0);
    for (std::size_t index = 0; index < profileCounts.size(); ++index) {
        prefix[index + 1] = prefix[index] + static_cast<double>(profileCounts[index]);
        result.profiledPulses += profileCounts[index];
    }
    if (result.profiledPulses == 0) {
        return result;
    }

    const int halfWindowBins = std::max(1, static_cast<int>(std::llround(0.5 * config.detrendWindowNs / config.profileBinWidthNs)));
    for (std::size_t index = 0; index < profileCounts.size(); ++index) {
        const std::size_t start = index > static_cast<std::size_t>(halfWindowBins)
                                ? index - static_cast<std::size_t>(halfWindowBins)
                                : 0;
        const std::size_t stop = std::min(profileCounts.size(), index + static_cast<std::size_t>(halfWindowBins) + 1);
        const double localMean = (prefix[stop] - prefix[start]) / static_cast<double>(stop - start);
        result.detrendedCounts[index] = static_cast<double>(profileCounts[index]) - localMean;
    }

    const int minLagBins = std::max(1, static_cast<int>(std::floor(config.periodSearchMinNs / config.profileBinWidthNs)));
    const int maxLagBins = std::max(minLagBins + 1, static_cast<int>(std::ceil(config.periodSearchMaxNs / config.profileBinWidthNs)));
    std::vector<PeriodCandidate> candidates;

    for (int lagBins = minLagBins; lagBins <= maxLagBins; ++lagBins) {
        double numerator = 0.0;
        double leftNorm = 0.0;
        double rightNorm = 0.0;
        for (std::size_t index = 0; index + static_cast<std::size_t>(lagBins) < result.detrendedCounts.size(); ++index) {
            const double left = result.detrendedCounts[index];
            const double right = result.detrendedCounts[index + static_cast<std::size_t>(lagBins)];
            numerator += left * right;
            leftNorm += left * left;
            rightNorm += right * right;
        }
        if (leftNorm <= 0.0 || rightNorm <= 0.0) {
            continue;
        }
        candidates.push_back(PeriodCandidate{
            lagBins * config.profileBinWidthNs,
            numerator / std::sqrt(leftNorm * rightNorm),
            lagBins,
        });
    }

    if (candidates.empty()) {
        return result;
    }

    std::sort(candidates.begin(), candidates.end(), [](const PeriodCandidate& left, const PeriodCandidate& right) {
        return left.score > right.score;
    });

    const std::size_t keepCount = std::min<std::size_t>(10, candidates.size());
    result.topCandidates.assign(candidates.begin(), candidates.begin() + keepCount);
    result.bestPeriodNs = result.topCandidates.front().periodNs;
    result.bestScore = result.topCandidates.front().score;
    result.available = true;

    const int bestLagBins = result.topCandidates.front().lagBins;
    if (bestLagBins > minLagBins && bestLagBins < maxLagBins) {
        auto find_score = [&](int lagBins) {
            for (const auto& candidate : candidates) {
                if (candidate.lagBins == lagBins) {
                    return candidate.score;
                }
            }
            return 0.0;
        };
        const double scoreLeft = find_score(bestLagBins - 1);
        const double scoreCenter = find_score(bestLagBins);
        const double scoreRight = find_score(bestLagBins + 1);
        const double denominator = scoreLeft - 2.0 * scoreCenter + scoreRight;
        if (std::abs(denominator) > 1.0e-12) {
            const double offsetBins = 0.5 * (scoreLeft - scoreRight) / denominator;
            if (std::abs(offsetBins) <= 1.0) {
                result.bestPeriodNs = (static_cast<double>(bestLagBins) + offsetBins) * config.profileBinWidthNs;
            }
        }
    }

    return result;
}

void write_delta_profile_file(const Config& config,
                              const std::vector<std::uint64_t>& profileCounts,
                              const std::vector<double>& detrendedCounts) {
    std::ofstream output(config.deltaProfileOutputPath);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to open delta profile output file: " + config.deltaProfileOutputPath);
    }

    output << "# columns: delta_ns count detrended_count\n";
    output << "# profile_bin_width_ns=" << format_double(config.profileBinWidthNs) << '\n';
    output << "# profile_tot_gate_ns=[" << format_double(config.profileTotMinNs) << ", "
           << format_double(config.profileTotMaxNs) << ")\n";
    output << "# detrend_window_ns=" << format_double(config.detrendWindowNs) << '\n';
    for (std::size_t index = 0; index < profileCounts.size(); ++index) {
        const double deltaCenterNs = config.deltaMinNs + (static_cast<double>(index) + 0.5) * config.profileBinWidthNs;
        const double detrended = index < detrendedCounts.size() ? detrendedCounts[index] : 0.0;
        output << std::fixed << std::setprecision(6)
               << deltaCenterNs << ' ' << profileCounts[index] << ' ' << detrended << '\n';
    }
}

void write_period_summary_file(const Config& config,
                               const PeriodEstimate& estimate) {
    std::ofstream output(config.periodSummaryOutputPath);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to open period summary output file: " + config.periodSummaryOutputPath);
    }

    output << "delta_range_ns=[" << format_double(config.deltaMinNs) << ", "
           << format_double(config.deltaMaxNs) << ")\n";
    output << "max_parent_delta_ns=" << format_double(config.maxParentDeltaNs) << '\n';
    output << "profile_bin_width_ns=" << format_double(config.profileBinWidthNs) << '\n';
    output << "profile_tot_gate_ns=[" << format_double(config.profileTotMinNs) << ", "
           << format_double(config.profileTotMaxNs) << ")\n";
    output << "detrend_window_ns=" << format_double(config.detrendWindowNs) << '\n';
    output << "period_search_ns=[" << format_double(config.periodSearchMinNs) << ", "
           << format_double(config.periodSearchMaxNs) << "]\n";
    output << "profiled_pulses=" << estimate.profiledPulses << '\n';
    output << "period_estimate_available=" << (estimate.available ? "true" : "false") << '\n';
    if (!estimate.available) {
        return;
    }

    output << "best_period_ns=" << format_double(estimate.bestPeriodNs) << '\n';
    output << "best_period_score=" << format_double(estimate.bestScore) << "\n\n";
    output << "top_period_candidates\n";
    for (const auto& candidate : estimate.topCandidates) {
        output << "  period_ns=" << format_double(candidate.periodNs)
               << " score=" << format_double(candidate.score)
               << " lag_bins=" << candidate.lagBins << '\n';
    }
}

std::string make_event_text(const std::string& type,
                            const std::string& streamId,
                            double triggerSeconds,
                            int channel,
                            const std::string& detail) {
    std::ostringstream stream;
    stream << type << " stream=" << streamId
           << " trigger_s=" << std::fixed << std::setprecision(9) << triggerSeconds
           << " ch=" << channel << ' ' << detail;
    return stream.str();
}

void write_summary_file(const Config& config,
                        const AnalysisStats& stats,
                        double elapsedSeconds) {
    std::ofstream output(config.summaryOutputPath);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to open summary output file: " + config.summaryOutputPath);
    }

    output << "source_file=" << config.inputPath << '\n';
    output << "reference_channel=" << config.referenceChannel << '\n';
    output << "signal_channel=" << config.signalChannel << '\n';
    output << "absolute_time_formula_ns=GlobalTriggerTime_seconds*1e9 + leading_time_ns\n";
    output << "parent_matching=closest earlier resolved reference pulse within the same TDC stream\n";
    output << "max_parent_delta_ns=" << format_double(config.maxParentDeltaNs) << '\n';
    output << "x_axis=delta_to_parent_ch0_no_modulo_ns\n";
    output << "modulo_wrap=not_used\n";
    output << "histogram_output=" << config.histogramOutputPath << '\n';
    output << "delta_range_ns=[" << format_double(config.deltaMinNs) << ", "
           << format_double(config.deltaMaxNs) << ")\n";
    output << "tot_range_ns=[" << format_double(config.totMinNs) << ", "
           << format_double(config.totMaxNs) << ")\n";
    output << "delta_bins=" << config.deltaBins << '\n';
    output << "tot_bins=" << config.totBins << '\n';
    output << "completed_tdc_blocks=" << stats.completedBlocks << '\n';
    output << "elapsed_seconds=" << format_double(elapsedSeconds) << "\n\n";

    const auto write_channel = [&](const std::string& label, const PulseStats& channelStats) {
        output << label << "\n";
        output << "  rows=" << channelStats.rows << '\n';
        output << "  rises=" << channelStats.rises << '\n';
        output << "  falls=" << channelStats.falls << '\n';
        output << "  sentinel_rows=" << channelStats.sentinelRows << '\n';
        output << "  duplicate_rise=" << channelStats.duplicateRise << '\n';
        output << "  orphan_rise_end_block=" << channelStats.orphanRiseEndBlock << '\n';
        output << "  orphan_fall_without_rise=" << channelStats.orphanFallWithoutRise << '\n';
        output << "  negative_tot=" << channelStats.negativeTot << '\n';
        output << "  zero_tot=" << channelStats.zeroTot << '\n';
        output << "  valid_positive_pairs=" << channelStats.validPositivePairs << '\n';
        if (channelStats.validPositivePairs > 0) {
            output << "  valid_positive_tot_min_ns=" << format_double(channelStats.validTotMin) << '\n';
            output << "  valid_positive_tot_max_ns=" << format_double(channelStats.validTotMax) << '\n';
            output << "  valid_positive_tot_mean_ns="
                   << format_double(channelStats.validTotSum / static_cast<double>(channelStats.validPositivePairs)) << '\n';
        }
        output << '\n';
    };

    write_channel("reference_channel_stats", stats.reference);
    write_channel("signal_channel_stats", stats.signal);

    output << "accepted_signal_pulses=" << stats.acceptedSignalPulses << '\n';
    output << "histogrammed_signal_pulses=" << stats.histogrammedSignalPulses << '\n';
    output << "delta_overflow=" << stats.deltaOverflow << '\n';
    output << "histogram_tot_overflow=" << stats.histogramTotOverflow << '\n';
    output << "rejected_no_earlier_parent=" << stats.rejectedNoEarlierParent << '\n';
    output << "rejected_bad_parent=" << stats.rejectedBadParent << '\n';
    output << "rejected_parent_too_old=" << stats.rejectedParentTooOld << '\n';
    output << "rejected_negative_parent_delta=" << stats.rejectedNegativeParentDelta << '\n';
    if (stats.accepted.count > 0) {
        output << "accepted_delta_min_ns=" << format_double(stats.accepted.deltaMin) << '\n';
        output << "accepted_delta_max_ns=" << format_double(stats.accepted.deltaMax) << '\n';
        output << "accepted_delta_mean_ns="
               << format_double(stats.accepted.deltaSum / static_cast<double>(stats.accepted.count)) << '\n';
        output << "accepted_tot_min_ns=" << format_double(stats.accepted.totMin) << '\n';
        output << "accepted_tot_max_ns=" << format_double(stats.accepted.totMax) << '\n';
        output << "accepted_tot_mean_ns="
               << format_double(stats.accepted.totSum / static_cast<double>(stats.accepted.count)) << '\n';
    }
    output << '\n';

    output << "accepted_by_stream\n";
    for (const auto& [streamId, count] : stats.acceptedByStream) {
        output << "  " << streamId << '=' << count << '\n';
    }
    output << '\n';

    output << "rejected_by_reason\n";
    for (const auto& [reason, count] : stats.rejectedByReason) {
        output << "  " << reason << '=' << count << '\n';
    }
    output << '\n';

    output << "sample_reference_anomalies\n";
    for (const auto& sample : stats.sampleReferenceAnomalies) {
        output << "  - " << sample << '\n';
    }
    output << '\n';

    output << "sample_signal_anomalies\n";
    for (const auto& sample : stats.sampleSignalAnomalies) {
        output << "  - " << sample << '\n';
    }
    output << '\n';

    output << "sample_parent_rejections\n";
    for (const auto& sample : stats.sampleParentRejects) {
        output << "  - " << sample << '\n';
    }
}

void write_histogram_grid(const Config& config,
                          const std::vector<std::uint64_t>& counts) {
    std::ofstream output(config.histogramOutputPath);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to open histogram output file: " + config.histogramOutputPath);
    }

    const double deltaBinWidth = (config.deltaMaxNs - config.deltaMinNs) / static_cast<double>(config.deltaBins);
    const double totBinWidth = (config.totMaxNs - config.totMinNs) / static_cast<double>(config.totBins);

    output << "# columns: delta_to_parent_ch0_us ch2_tot_ns count\n";
    output << "# x_axis=delta_to_parent_ch0_no_modulo_us\n";
    output << "# delta_min_ns=" << format_double(config.deltaMinNs) << '\n';
    output << "# delta_max_ns=" << format_double(config.deltaMaxNs) << '\n';
    output << "# delta_bins=" << config.deltaBins << '\n';
    output << "# tot_min_ns=" << format_double(config.totMinNs) << '\n';
    output << "# tot_max_ns=" << format_double(config.totMaxNs) << '\n';
    output << "# tot_bins=" << config.totBins << '\n';

    for (int totIndex = 0; totIndex < config.totBins; ++totIndex) {
        const double totCenter = config.totMinNs + (static_cast<double>(totIndex) + 0.5) * totBinWidth;
        const auto rowOffset = static_cast<std::size_t>(totIndex) * static_cast<std::size_t>(config.deltaBins);
        for (int deltaIndex = 0; deltaIndex < config.deltaBins; ++deltaIndex) {
            const double deltaCenterNs = config.deltaMinNs + (static_cast<double>(deltaIndex) + 0.5) * deltaBinWidth;
            output << std::fixed << std::setprecision(6)
                   << (deltaCenterNs / 1000.0) << ' ' << totCenter << ' '
                   << counts[rowOffset + static_cast<std::size_t>(deltaIndex)] << '\n';
        }
        output << '\n';
    }
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
        std::unordered_map<std::string, StreamState> streams;
        std::vector<std::uint64_t> histogram(static_cast<std::size_t>(config.deltaBins) * static_cast<std::size_t>(config.totBins), 0);
        std::vector<std::uint64_t> deltaProfile;
        const double deltaBinWidth = (config.deltaMaxNs - config.deltaMinNs) / static_cast<double>(config.deltaBins);
        const double totBinWidth = (config.totMaxNs - config.totMinNs) / static_cast<double>(config.totBins);
        if (!config.deltaProfileOutputPath.empty() || !config.periodSummaryOutputPath.empty()) {
            const auto profileBins = static_cast<std::size_t>(std::ceil((config.deltaMaxNs - config.deltaMinNs) / config.profileBinWidthNs));
            deltaProfile.assign(profileBins, 0);
        }

        const auto accept_signal = [&](const std::string& streamId,
                                       const SignalPulse& pulse,
                                       const ResolvedParent& parent) {
            const double deltaNs = pulse.absoluteLeadingNs - parent.absoluteLeadingNs;
            if (deltaNs < 0.0) {
                stats.rejectedNegativeParentDelta += 1;
                stats.rejectedByReason["negative_parent_delta"] += 1;
                add_sample(stats.sampleParentRejects,
                           make_event_text("negative_parent_delta", streamId, pulse.triggerSeconds,
                                           config.signalChannel,
                                           "signal_leading_ns=" + format_double(pulse.leadingNs) +
                                           " parent_leading_ns=" + format_double(parent.leadingNs)));
                return;
            }

            if (deltaNs >= config.maxParentDeltaNs) {
                stats.rejectedParentTooOld += 1;
                stats.rejectedByReason["parent_too_old"] += 1;
                add_sample(stats.sampleParentRejects,
                           make_event_text("parent_too_old", streamId, pulse.triggerSeconds,
                                           config.signalChannel,
                                           "signal_leading_ns=" + format_double(pulse.leadingNs) +
                                           " signal_tot_ns=" + format_double(pulse.totNs) +
                                           " delta_ns=" + format_double(deltaNs) +
                                           " max_parent_delta_ns=" + format_double(config.maxParentDeltaNs)));
                return;
            }

            stats.acceptedSignalPulses += 1;
            stats.accepted.count += 1;
            stats.accepted.deltaMin = std::min(stats.accepted.deltaMin, deltaNs);
            stats.accepted.deltaMax = std::max(stats.accepted.deltaMax, deltaNs);
            stats.accepted.deltaSum += deltaNs;
            stats.accepted.totMin = std::min(stats.accepted.totMin, pulse.totNs);
            stats.accepted.totMax = std::max(stats.accepted.totMax, pulse.totNs);
            stats.accepted.totSum += pulse.totNs;
            stats.acceptedByStream[streamId] += 1;

            if (!deltaProfile.empty() &&
                pulse.totNs >= config.profileTotMinNs && pulse.totNs < config.profileTotMaxNs &&
                deltaNs >= config.deltaMinNs && deltaNs < config.deltaMaxNs) {
                const int profileIndex = std::min(static_cast<int>(deltaProfile.size()) - 1,
                                                  static_cast<int>(std::floor((deltaNs - config.deltaMinNs) / config.profileBinWidthNs)));
                deltaProfile[static_cast<std::size_t>(profileIndex)] += 1;
            }

            if (deltaNs < config.deltaMinNs || deltaNs >= config.deltaMaxNs) {
                stats.deltaOverflow += 1;
                return;
            }
            if (pulse.totNs < config.totMinNs || pulse.totNs >= config.totMaxNs) {
                stats.histogramTotOverflow += 1;
                return;
            }

            const int deltaIndex = std::min(config.deltaBins - 1,
                                            static_cast<int>(std::floor((deltaNs - config.deltaMinNs) / deltaBinWidth)));
            const int totIndex = std::min(config.totBins - 1,
                                          static_cast<int>(std::floor((pulse.totNs - config.totMinNs) / totBinWidth)));
            const auto flatIndex = static_cast<std::size_t>(totIndex) * static_cast<std::size_t>(config.deltaBins)
                                 + static_cast<std::size_t>(deltaIndex);
            histogram[flatIndex] += 1;
            stats.histogrammedSignalPulses += 1;
        };

        const auto reject_signal = [&](const std::string& reason,
                                       const std::string& streamId,
                                       const SignalPulse& pulse,
                                       const std::string& detail) {
            if (reason == "no_earlier_parent") {
                stats.rejectedNoEarlierParent += 1;
            } else {
                stats.rejectedBadParent += 1;
            }
            stats.rejectedByReason[reason] += 1;
            add_sample(stats.sampleParentRejects,
                       make_event_text(reason, streamId, pulse.triggerSeconds, config.signalChannel,
                                       "signal_leading_ns=" + format_double(pulse.leadingNs) +
                                       " signal_tot_ns=" + format_double(pulse.totNs) + ' ' + detail));
        };

        const auto find_reference_parent = [](const StreamState& stream,
                                              double signalAbsoluteLeadingNs) -> const ResolvedParent* {
            for (auto iterator = stream.resolvedReferenceHistory.rbegin();
                 iterator != stream.resolvedReferenceHistory.rend();
                 ++iterator) {
                if (iterator->absoluteLeadingNs <= signalAbsoluteLeadingNs) {
                    return &(*iterator);
                }
            }
            return nullptr;
        };

        const auto resolve_reference = [&](const std::string& streamId,
                                           StreamState& stream,
                                           bool usable,
                                           const std::string& reason,
                                           double totNs) {
            ResolvedParent resolved;
            resolved.exists = true;
            resolved.usable = usable;
            resolved.triggerSeconds = stream.pendingReference.triggerSeconds;
            resolved.leadingNs = stream.pendingReference.leadingNs;
            resolved.absoluteLeadingNs = stream.pendingReference.absoluteLeadingNs;
            resolved.totNs = totNs;
            resolved.firstCol = stream.pendingReference.firstCol;
            resolved.reason = reason;

            stream.resolvedReferenceHistory.push_back(resolved);
            if (stream.resolvedReferenceHistory.size() > kResolvedHistoryLimit) {
                stream.resolvedReferenceHistory.erase(stream.resolvedReferenceHistory.begin());
            }

            for (const auto& pulse : stream.waitingSignalPulses) {
                if (usable) {
                    accept_signal(streamId, pulse, resolved);
                } else {
                    reject_signal("bad_parent_" + reason, streamId, pulse,
                                  "parent_leading_ns=" + format_double(stream.pendingReference.leadingNs));
                }
            }
            stream.waitingSignalPulses.clear();
            stream.pendingReference = PendingPulse{};
        };

        const auto flush_block = [&](const std::string& streamId) {
            auto streamIt = streams.find(streamId);
            if (streamIt == streams.end()) {
                return;
            }
            auto& stream = streamIt->second;

            if (stream.pendingReference.hasValue) {
                stats.reference.orphanRiseEndBlock += 1;
                add_sample(stats.sampleReferenceAnomalies,
                           make_event_text("reference_block_end_orphan_rise", streamId,
                                           stream.pendingReference.triggerSeconds, config.referenceChannel,
                                           "leading_ns=" + format_double(stream.pendingReference.leadingNs)));
                resolve_reference(streamId, stream, false, "block_end_orphan_rise", 0.0);
            }

            if (stream.pendingSignal.hasValue) {
                stats.signal.orphanRiseEndBlock += 1;
                add_sample(stats.sampleSignalAnomalies,
                           make_event_text("signal_block_end_orphan_rise", streamId,
                                           stream.pendingSignal.triggerSeconds, config.signalChannel,
                                           "leading_ns=" + format_double(stream.pendingSignal.leadingNs)));
                stream.pendingSignal = PendingPulse{};
            }
        };

        std::string currentStreamId = "UNKNOWN";
        std::optional<double> currentTrigger;
        std::string line;
        const auto startedAt = std::chrono::steady_clock::now();

        while (std::getline(input, line)) {
            if (starts_with(line, "TDC ")) {
                if (currentStreamId != "UNKNOWN") {
                    flush_block(currentStreamId);
                    stats.completedBlocks += 1;
                }
                currentStreamId = parse_tdc_id(line);
                currentTrigger.reset();
                continue;
            }

            if (const auto trigger = parse_trigger(line); trigger.has_value()) {
                currentTrigger = trigger;
                continue;
            }

            int firstCol = 0;
            int channel = 0;
            int edge = 0;
            double valueNs = 0.0;
            if (!parse_row(line, firstCol, channel, edge, valueNs)) {
                continue;
            }
            if (!currentTrigger.has_value()) {
                continue;
            }

            if (channel != config.referenceChannel && channel != config.signalChannel) {
                continue;
            }

            auto& stream = streams[currentStreamId];
            PulseStats& channelStats = channel == config.referenceChannel ? stats.reference : stats.signal;
            channelStats.rows += 1;

            if (valueNs == kSentinelValue) {
                channelStats.sentinelRows += 1;
                add_sample(channel == config.referenceChannel ? stats.sampleReferenceAnomalies : stats.sampleSignalAnomalies,
                           make_event_text(channel == config.referenceChannel ? "reference_sentinel_row" : "signal_sentinel_row",
                                           currentStreamId, *currentTrigger, channel,
                                           std::string("edge=") + std::to_string(edge)));
                continue;
            }

            PendingPulse& pending = channel == config.referenceChannel ? stream.pendingReference : stream.pendingSignal;
            if (edge == 1) {
                channelStats.rises += 1;
                if (pending.hasValue) {
                    channelStats.duplicateRise += 1;
                    add_sample(channel == config.referenceChannel ? stats.sampleReferenceAnomalies : stats.sampleSignalAnomalies,
                               make_event_text(channel == config.referenceChannel ? "reference_duplicate_rise" : "signal_duplicate_rise",
                                               currentStreamId, *currentTrigger, channel,
                                               "previous_leading_ns=" + format_double(pending.leadingNs) +
                                               " new_leading_ns=" + format_double(valueNs)));
                    if (channel == config.referenceChannel && !stream.waitingSignalPulses.empty()) {
                        for (const auto& pulse : stream.waitingSignalPulses) {
                            reject_signal("bad_parent_duplicate_rise_before_fall", currentStreamId, pulse,
                                          "parent_leading_ns=" + format_double(pending.leadingNs));
                        }
                        stream.waitingSignalPulses.clear();
                    }
                }
                pending.hasValue = true;
                pending.triggerSeconds = *currentTrigger;
                pending.leadingNs = valueNs;
                pending.absoluteLeadingNs = absolute_time_ns(*currentTrigger, valueNs);
                pending.firstCol = firstCol;
                continue;
            }

            channelStats.falls += 1;
            if (!pending.hasValue) {
                channelStats.orphanFallWithoutRise += 1;
                add_sample(channel == config.referenceChannel ? stats.sampleReferenceAnomalies : stats.sampleSignalAnomalies,
                           make_event_text(channel == config.referenceChannel ? "reference_fall_without_rise" : "signal_fall_without_rise",
                                           currentStreamId, *currentTrigger, channel,
                                           "fall_ns=" + format_double(valueNs)));
                continue;
            }

            const double totNs = valueNs - pending.leadingNs;
            if (totNs < 0.0) {
                channelStats.negativeTot += 1;
                add_sample(channel == config.referenceChannel ? stats.sampleReferenceAnomalies : stats.sampleSignalAnomalies,
                           make_event_text(channel == config.referenceChannel ? "reference_negative_tot" : "signal_negative_tot",
                                           currentStreamId, *currentTrigger, channel,
                                           "leading_ns=" + format_double(pending.leadingNs) +
                                           " fall_ns=" + format_double(valueNs) +
                                           " tot_ns=" + format_double(totNs)));
                if (channel == config.referenceChannel) {
                    resolve_reference(currentStreamId, stream, false, "negative_tot", totNs);
                } else {
                    pending = PendingPulse{};
                }
                continue;
            }

            if (totNs == 0.0) {
                channelStats.zeroTot += 1;
                add_sample(channel == config.referenceChannel ? stats.sampleReferenceAnomalies : stats.sampleSignalAnomalies,
                           make_event_text(channel == config.referenceChannel ? "reference_zero_tot" : "signal_zero_tot",
                                           currentStreamId, *currentTrigger, channel,
                                           "leading_ns=" + format_double(pending.leadingNs) +
                                           " fall_ns=" + format_double(valueNs)));
                if (channel == config.referenceChannel) {
                    resolve_reference(currentStreamId, stream, false, "zero_tot", totNs);
                } else {
                    pending = PendingPulse{};
                }
                continue;
            }

            update_positive_tot(channelStats, totNs);
            if (channel == config.referenceChannel) {
                resolve_reference(currentStreamId, stream, true, "ok", totNs);
                continue;
            }

            const SignalPulse pulse{*currentTrigger, pending.leadingNs, pending.absoluteLeadingNs, totNs, pending.firstCol};
            pending = PendingPulse{};

            if (stream.pendingReference.hasValue && pulse.absoluteLeadingNs >= stream.pendingReference.absoluteLeadingNs) {
                stream.waitingSignalPulses.push_back(pulse);
                continue;
            }

            const ResolvedParent* parent = find_reference_parent(stream, pulse.absoluteLeadingNs);
            if (parent == nullptr) {
                reject_signal("no_earlier_parent", currentStreamId, pulse, "reason=no_reference_seen_yet");
                continue;
            }

            if (!parent->usable) {
                reject_signal("bad_parent_" + parent->reason,
                              currentStreamId, pulse,
                              "parent_leading_ns=" + format_double(parent->leadingNs));
                continue;
            }

            accept_signal(currentStreamId, pulse, *parent);
        }

        if (currentStreamId != "UNKNOWN") {
            flush_block(currentStreamId);
            stats.completedBlocks += 1;
        }

        const auto endedAt = std::chrono::steady_clock::now();
        const double elapsedSeconds = std::chrono::duration_cast<std::chrono::duration<double>>(endedAt - startedAt).count();

        PeriodEstimate periodEstimate;
        if (!deltaProfile.empty()) {
            periodEstimate = estimate_period(config, deltaProfile);
            if (!config.deltaProfileOutputPath.empty()) {
                write_delta_profile_file(config, deltaProfile, periodEstimate.detrendedCounts);
            }
            if (!config.periodSummaryOutputPath.empty()) {
                write_period_summary_file(config, periodEstimate);
            }
        }

        write_histogram_grid(config, histogram);
        write_summary_file(config, stats, elapsedSeconds);

        std::cout << "Source file: " << config.inputPath << '\n';
        std::cout << "Summary output: " << config.summaryOutputPath << '\n';
        std::cout << "Histogram grid: " << config.histogramOutputPath << '\n';
        if (!config.deltaProfileOutputPath.empty()) {
            std::cout << "Delta profile: " << config.deltaProfileOutputPath << '\n';
        }
        if (!config.periodSummaryOutputPath.empty()) {
            std::cout << "Period summary: " << config.periodSummaryOutputPath << '\n';
            if (periodEstimate.available) {
                std::cout << "Estimated period [ns]: " << format_double(periodEstimate.bestPeriodNs) << '\n';
                std::cout << "Estimated period score: " << format_double(periodEstimate.bestScore) << '\n';
            }
        }
        std::cout << "Completed TDC blocks: " << stats.completedBlocks << '\n';
        std::cout << "Elapsed seconds: " << format_double(elapsedSeconds) << '\n';
        std::cout << "Accepted signal pulses: " << stats.acceptedSignalPulses << '\n';
        std::cout << "Histogrammed signal pulses: " << stats.histogrammedSignalPulses << '\n';
        std::cout << "Delta overflow: " << stats.deltaOverflow << '\n';
        std::cout << "Histogram TOT overflow: " << stats.histogramTotOverflow << '\n';
        std::cout << "Rejected no earlier parent: " << stats.rejectedNoEarlierParent << '\n';
        std::cout << "Rejected bad parent: " << stats.rejectedBadParent << '\n';
        std::cout << "Rejected negative parent delta: " << stats.rejectedNegativeParentDelta << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}