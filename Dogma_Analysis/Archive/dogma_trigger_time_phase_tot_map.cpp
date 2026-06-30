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

struct Config {
    std::string inputPath;
    std::string summaryOutputPath;
    std::string histogramOutputPath;
    int signalChannel = 2;
    double phasePeriodNs = 40.0;
    int phaseBins = 400;
    double totMinNs = 0.0;
    double totMaxNs = 128.0;
    int totBins = 512;
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
    double phaseMin = std::numeric_limits<double>::infinity();
    double phaseMax = -std::numeric_limits<double>::infinity();
    double phaseSum = 0.0;
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

struct StreamState {
    PendingPulse pendingSignal;
};

struct AnalysisStats {
    PulseStats signal;
    AcceptedStats accepted;
    std::uint64_t completedBlocks = 0;
    std::uint64_t histogrammedSignalPulses = 0;
    std::uint64_t histogramTotOverflow = 0;
    std::map<std::string, std::uint64_t> acceptedByStream;
    std::vector<std::string> sampleSignalAnomalies;
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
        } else if (arg == "--signal-channel" && index + 1 < argc) {
            config.signalChannel = std::stoi(argv[++index]);
        } else if (arg == "--phase-period-ns" && index + 1 < argc) {
            config.phasePeriodNs = std::stod(argv[++index]);
        } else if (arg == "--phase-bins" && index + 1 < argc) {
            config.phaseBins = std::stoi(argv[++index]);
        } else if (arg == "--tot-min-ns" && index + 1 < argc) {
            config.totMinNs = std::stod(argv[++index]);
        } else if (arg == "--tot-max-ns" && index + 1 < argc) {
            config.totMaxNs = std::stod(argv[++index]);
        } else if (arg == "--tot-bins" && index + 1 < argc) {
            config.totBins = std::stoi(argv[++index]);
        } else if (arg == "--help") {
            std::cout
                << "Usage: dogma_trigger_time_phase_tot_map --input <file> --summary-output <txt> --histogram-output <tsv>\n"
                << "  [--signal-channel 2] [--phase-period-ns 40.0] [--phase-bins 400]\n"
                << "  [--tot-min-ns 0] [--tot-max-ns 128] [--tot-bins 512]\n";
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
    if (config.phaseBins <= 0 || config.totBins <= 0) {
        throw std::runtime_error("Histogram bin counts must be positive");
    }
    if (config.phasePeriodNs <= 0.0) {
        throw std::runtime_error("Phase period must be positive");
    }
    if (config.totMaxNs <= config.totMinNs) {
        throw std::runtime_error("tot-max-ns must be larger than tot-min-ns");
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

double positive_mod(double value, double period) {
    const double wrapped = std::fmod(value, period);
    return wrapped < 0.0 ? wrapped + period : wrapped;
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
    output << "signal_channel=" << config.signalChannel << '\n';
    output << "phase_reference=trigger_time_referenced\n";
    output << "absolute_time_formula_ns=GlobalTriggerTime_seconds*1e9 + leading_time_ns\n";
    output << "phase_formula_ns=(signal_leading_abs_ns) mod " << format_double(config.phasePeriodNs) << '\n';
    output << "histogram_output=" << config.histogramOutputPath << '\n';
    output << "phase_bins=" << config.phaseBins << '\n';
    output << "tot_bins=" << config.totBins << '\n';
    output << "tot_range_ns=[" << format_double(config.totMinNs) << ", "
           << format_double(config.totMaxNs) << ")\n";
    output << "completed_tdc_blocks=" << stats.completedBlocks << '\n';
    output << "elapsed_seconds=" << format_double(elapsedSeconds) << "\n\n";

    output << "signal_channel_stats\n";
    output << "  rows=" << stats.signal.rows << '\n';
    output << "  rises=" << stats.signal.rises << '\n';
    output << "  falls=" << stats.signal.falls << '\n';
    output << "  sentinel_rows=" << stats.signal.sentinelRows << '\n';
    output << "  duplicate_rise=" << stats.signal.duplicateRise << '\n';
    output << "  orphan_rise_end_block=" << stats.signal.orphanRiseEndBlock << '\n';
    output << "  orphan_fall_without_rise=" << stats.signal.orphanFallWithoutRise << '\n';
    output << "  negative_tot=" << stats.signal.negativeTot << '\n';
    output << "  zero_tot=" << stats.signal.zeroTot << '\n';
    output << "  valid_positive_pairs=" << stats.signal.validPositivePairs << '\n';
    if (stats.signal.validPositivePairs > 0) {
        output << "  valid_positive_tot_min_ns=" << format_double(stats.signal.validTotMin) << '\n';
        output << "  valid_positive_tot_max_ns=" << format_double(stats.signal.validTotMax) << '\n';
        output << "  valid_positive_tot_mean_ns="
               << format_double(stats.signal.validTotSum / static_cast<double>(stats.signal.validPositivePairs)) << '\n';
    }
    output << '\n';

    output << "accepted_signal_pulses=" << stats.accepted.count << '\n';
    output << "histogrammed_signal_pulses=" << stats.histogrammedSignalPulses << '\n';
    output << "histogram_tot_overflow=" << stats.histogramTotOverflow << '\n';
    if (stats.accepted.count > 0) {
        output << "accepted_phase_min_ns=" << format_double(stats.accepted.phaseMin) << '\n';
        output << "accepted_phase_max_ns=" << format_double(stats.accepted.phaseMax) << '\n';
        output << "accepted_phase_mean_ns="
               << format_double(stats.accepted.phaseSum / static_cast<double>(stats.accepted.count)) << '\n';
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

    output << "sample_signal_anomalies\n";
    for (const auto& sample : stats.sampleSignalAnomalies) {
        output << "  - " << sample << '\n';
    }
}

void write_histogram_file(const Config& config,
                          const std::vector<std::uint64_t>& counts) {
    std::ofstream output(config.histogramOutputPath);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to open histogram output file: " + config.histogramOutputPath);
    }

    output << "# phase_bins=" << config.phaseBins << '\n';
    output << "# phase_min_ns=0\n";
    output << "# phase_max_ns=" << format_double(config.phasePeriodNs) << '\n';
    output << "# tot_bins=" << config.totBins << '\n';
    output << "# tot_min_ns=" << format_double(config.totMinNs) << '\n';
    output << "# tot_max_ns=" << format_double(config.totMaxNs) << '\n';
    for (int totIndex = 0; totIndex < config.totBins; ++totIndex) {
        const auto rowOffset = static_cast<std::size_t>(totIndex) * static_cast<std::size_t>(config.phaseBins);
        for (int phaseIndex = 0; phaseIndex < config.phaseBins; ++phaseIndex) {
            if (phaseIndex > 0) {
                output << '\t';
            }
            output << counts[rowOffset + static_cast<std::size_t>(phaseIndex)];
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
        std::vector<std::uint64_t> histogram(static_cast<std::size_t>(config.phaseBins) * static_cast<std::size_t>(config.totBins), 0);
        const double phaseBinWidth = config.phasePeriodNs / static_cast<double>(config.phaseBins);
        const double totBinWidth = (config.totMaxNs - config.totMinNs) / static_cast<double>(config.totBins);

        const auto accept_signal = [&](const std::string& streamId,
                                       const PendingPulse& pulse,
                                       double totNs) {
            const double phaseNs = positive_mod(pulse.absoluteLeadingNs, config.phasePeriodNs);
            stats.accepted.count += 1;
            stats.accepted.phaseMin = std::min(stats.accepted.phaseMin, phaseNs);
            stats.accepted.phaseMax = std::max(stats.accepted.phaseMax, phaseNs);
            stats.accepted.phaseSum += phaseNs;
            stats.accepted.totMin = std::min(stats.accepted.totMin, totNs);
            stats.accepted.totMax = std::max(stats.accepted.totMax, totNs);
            stats.accepted.totSum += totNs;
            stats.acceptedByStream[streamId] += 1;

            if (totNs < config.totMinNs || totNs >= config.totMaxNs) {
                stats.histogramTotOverflow += 1;
                return;
            }

            const int phaseIndex = std::min(config.phaseBins - 1,
                                            static_cast<int>(std::floor(phaseNs / phaseBinWidth)));
            const int totIndex = std::min(config.totBins - 1,
                                          static_cast<int>(std::floor((totNs - config.totMinNs) / totBinWidth)));
            const auto flatIndex = static_cast<std::size_t>(totIndex) * static_cast<std::size_t>(config.phaseBins)
                                 + static_cast<std::size_t>(phaseIndex);
            histogram[flatIndex] += 1;
            stats.histogrammedSignalPulses += 1;
        };

        const auto flush_block = [&](const std::string& streamId) {
            auto streamIt = streams.find(streamId);
            if (streamIt == streams.end()) {
                return;
            }
            auto& pending = streamIt->second.pendingSignal;
            if (!pending.hasValue) {
                return;
            }
            stats.signal.orphanRiseEndBlock += 1;
            add_sample(stats.sampleSignalAnomalies,
                       make_event_text("signal_block_end_orphan_rise", streamId,
                                       pending.triggerSeconds, config.signalChannel,
                                       "leading_ns=" + format_double(pending.leadingNs)));
            pending = PendingPulse{};
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
            if (!currentTrigger.has_value() || channel != config.signalChannel) {
                continue;
            }

            auto& pending = streams[currentStreamId].pendingSignal;
            stats.signal.rows += 1;

            if (valueNs == kSentinelValue) {
                stats.signal.sentinelRows += 1;
                add_sample(stats.sampleSignalAnomalies,
                           make_event_text("signal_sentinel_row", currentStreamId,
                                           *currentTrigger, channel,
                                           std::string("edge=") + std::to_string(edge)));
                continue;
            }

            if (edge == 1) {
                stats.signal.rises += 1;
                if (pending.hasValue) {
                    stats.signal.duplicateRise += 1;
                    add_sample(stats.sampleSignalAnomalies,
                               make_event_text("signal_duplicate_rise", currentStreamId,
                                               *currentTrigger, channel,
                                               "previous_leading_ns=" + format_double(pending.leadingNs) +
                                               " new_leading_ns=" + format_double(valueNs)));
                }
                pending.hasValue = true;
                pending.triggerSeconds = *currentTrigger;
                pending.leadingNs = valueNs;
                pending.absoluteLeadingNs = absolute_time_ns(*currentTrigger, valueNs);
                pending.firstCol = firstCol;
                continue;
            }

            stats.signal.falls += 1;
            if (!pending.hasValue) {
                stats.signal.orphanFallWithoutRise += 1;
                add_sample(stats.sampleSignalAnomalies,
                           make_event_text("signal_fall_without_rise", currentStreamId,
                                           *currentTrigger, channel,
                                           "fall_ns=" + format_double(valueNs)));
                continue;
            }

            const double totNs = valueNs - pending.leadingNs;
            if (totNs < 0.0) {
                stats.signal.negativeTot += 1;
                add_sample(stats.sampleSignalAnomalies,
                           make_event_text("signal_negative_tot", currentStreamId,
                                           *currentTrigger, channel,
                                           "leading_ns=" + format_double(pending.leadingNs) +
                                           " fall_ns=" + format_double(valueNs) +
                                           " tot_ns=" + format_double(totNs)));
                pending = PendingPulse{};
                continue;
            }

            if (totNs == 0.0) {
                stats.signal.zeroTot += 1;
                add_sample(stats.sampleSignalAnomalies,
                           make_event_text("signal_zero_tot", currentStreamId,
                                           *currentTrigger, channel,
                                           "leading_ns=" + format_double(pending.leadingNs) +
                                           " fall_ns=" + format_double(valueNs)));
                pending = PendingPulse{};
                continue;
            }

            update_positive_tot(stats.signal, totNs);
            accept_signal(currentStreamId, pending, totNs);
            pending = PendingPulse{};
        }

        if (currentStreamId != "UNKNOWN") {
            flush_block(currentStreamId);
            stats.completedBlocks += 1;
        }

        const auto endedAt = std::chrono::steady_clock::now();
        const double elapsedSeconds = std::chrono::duration_cast<std::chrono::duration<double>>(endedAt - startedAt).count();

        write_histogram_file(config, histogram);
        write_summary_file(config, stats, elapsedSeconds);

        std::cout << "Source file: " << config.inputPath << '\n';
        std::cout << "Summary output: " << config.summaryOutputPath << '\n';
        std::cout << "Histogram TSV: " << config.histogramOutputPath << '\n';
        std::cout << "Completed TDC blocks: " << stats.completedBlocks << '\n';
        std::cout << "Elapsed seconds: " << format_double(elapsedSeconds) << '\n';
        std::cout << "Accepted signal pulses: " << stats.accepted.count << '\n';
        std::cout << "Histogrammed signal pulses: " << stats.histogrammedSignalPulses << '\n';
        std::cout << "Histogram TOT overflow: " << stats.histogramTotOverflow << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}