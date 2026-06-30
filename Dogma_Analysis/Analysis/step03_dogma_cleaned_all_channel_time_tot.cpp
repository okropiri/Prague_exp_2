#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr int kChannelCount = 32;

struct Config {
    std::string inputPath;
    std::string outputRoot;
    std::string runKey;
    std::string triggerRefDirName = "Trigger_ref_ToT";
    std::string ch0RefDirName = "Ch0_ref_TOT";
    std::string totDistribDirName = "TOT_distrib";
    double triggerRefMinNs = -6000.0;
    double triggerRefMaxNs = 6000.0;
    int triggerRefBins = 2400;
    double ch0RefMinNs = -6000.0;
    double ch0RefMaxNs = 6000.0;
    int ch0RefBins = 2400;
    double totMinNs = 0.0;
    double totMaxNs = 128.0;
    int totBins = 1280;
    double ch0ValidRiseMinNs = -410.0;
    double ch0ValidRiseMaxNs = -395.0;
    double ch0ValidTotMinNs = 16.5;
    double ch0ValidTotMaxNs = 19.5;
};

struct PulseRow {
    std::uint64_t windowIndex = 0;
    int channel = -1;
    double riseNs = 0.0;
    double totNs = 0.0;
};

struct PulseHit {
    double riseNs = 0.0;
    double totNs = 0.0;
};

struct WindowState {
    bool hasValue = false;
    std::uint64_t windowIndex = 0;
    std::array<std::vector<PulseHit>, kChannelCount> hitsByChannel;
};

struct Stats {
    std::string sourceInputFile;
    std::uint64_t parsedPulseRows = 0;
    std::uint64_t invalidRows = 0;
    std::uint64_t windowsSeen = 0;
    std::uint64_t validCh0Windows = 0;
    std::uint64_t windowsWithoutCh0 = 0;
    std::uint64_t windowsWithMultipleCh0 = 0;
    std::uint64_t windowsWithoutValidCh0Candidate = 0;
    std::uint64_t windowsWithExactlyOneValidCh0Candidate = 0;
    std::uint64_t windowsWithMultipleValidCh0Candidates = 0;
    std::uint64_t totalCh0Candidates = 0;
    std::uint64_t totalValidCh0Candidates = 0;
    std::uint64_t rejectedCh0OutsideTimeRange = 0;
    std::uint64_t rejectedCh0OutsideTotRange = 0;
    std::uint64_t triggerRefAcceptedPulses = 0;
    std::uint64_t triggerRefRejectedOutsideTimeRange = 0;
    std::uint64_t triggerRefRejectedOutsideTotRange = 0;
    std::uint64_t triggerRefRejectedOutsideTimeAndTotRange = 0;
    std::uint64_t ch0RefAcceptedPulses = 0;
    std::uint64_t ch0RefRejectedOutsideTimeRange = 0;
    std::uint64_t ch0RefRejectedOutsideTotRange = 0;
    std::uint64_t ch0RefRejectedOutsideTimeAndTotRange = 0;
    double observedTotMinNs = std::numeric_limits<double>::max();
    double observedTotMaxNs = std::numeric_limits<double>::lowest();
    std::uint64_t wholeRunTotAcceptedPulses = 0;
    std::uint64_t wholeRunTotRejectedBelowMin = 0;
    std::array<std::uint64_t, kChannelCount> pulseCountByChannel{};
    std::array<std::uint64_t, kChannelCount> triggerRefAcceptedByChannel{};
    std::array<std::uint64_t, kChannelCount> ch0RefAcceptedByChannel{};
    std::array<std::uint64_t, kChannelCount> wholeRunTotAcceptedByChannel{};
};

using SparseHistogram = std::array<std::unordered_map<std::uint64_t, std::uint64_t>, kChannelCount>;
using TotDistribution = std::array<std::vector<std::uint64_t>, kChannelCount>;

bool starts_with(std::string_view text, std::string_view prefix) {
    return text.size() >= prefix.size() && text.substr(0, prefix.size()) == prefix;
}

std::string trim(std::string_view text) {
    std::size_t start = 0;
    while (start < text.size() && std::isspace(static_cast<unsigned char>(text[start]))) {
        start += 1;
    }
    std::size_t end = text.size();
    while (end > start && std::isspace(static_cast<unsigned char>(text[end - 1]))) {
        end -= 1;
    }
    return std::string(text.substr(start, end - start));
}

std::string format_double(double value) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(6) << value;
    std::string text = stream.str();
    while (!text.empty() && text.back() == '0') {
        text.pop_back();
    }
    if (!text.empty() && text.back() == '.') {
        text.pop_back();
    }
    return text.empty() ? "0" : text;
}

std::string strip_known_suffixes(std::string text) {
    bool changed = true;
    while (changed) {
        changed = false;
        for (const std::string_view suffix : {std::string_view{"_rawRefined_pulses.tsv"}, std::string_view{".tsv"}, std::string_view{".txt"}}) {
            if (text.size() >= suffix.size() && text.substr(text.size() - suffix.size()) == suffix) {
                text.erase(text.size() - suffix.size());
                changed = true;
                break;
            }
        }
    }
    return text;
}

Config parse_args(int argc, char** argv) {
    Config config;
    for (int index = 1; index < argc; ++index) {
        const std::string arg = argv[index];
        if (arg == "--input" && index + 1 < argc) {
            config.inputPath = argv[++index];
        } else if (arg == "--output-root" && index + 1 < argc) {
            config.outputRoot = argv[++index];
        } else if (arg == "--run-key" && index + 1 < argc) {
            config.runKey = argv[++index];
        } else if (arg == "--trigger-ref-dir-name" && index + 1 < argc) {
            config.triggerRefDirName = argv[++index];
        } else if (arg == "--ch0ref-dir-name" && index + 1 < argc) {
            config.ch0RefDirName = argv[++index];
        } else if (arg == "--trigger-ref-min-ns" && index + 1 < argc) {
            config.triggerRefMinNs = std::stod(argv[++index]);
        } else if (arg == "--trigger-ref-max-ns" && index + 1 < argc) {
            config.triggerRefMaxNs = std::stod(argv[++index]);
        } else if (arg == "--trigger-ref-bins" && index + 1 < argc) {
            config.triggerRefBins = std::stoi(argv[++index]);
        } else if (arg == "--ch0ref-min-ns" && index + 1 < argc) {
            config.ch0RefMinNs = std::stod(argv[++index]);
        } else if (arg == "--ch0ref-max-ns" && index + 1 < argc) {
            config.ch0RefMaxNs = std::stod(argv[++index]);
        } else if (arg == "--ch0ref-bins" && index + 1 < argc) {
            config.ch0RefBins = std::stoi(argv[++index]);
        } else if (arg == "--tot-min-ns" && index + 1 < argc) {
            config.totMinNs = std::stod(argv[++index]);
        } else if (arg == "--tot-max-ns" && index + 1 < argc) {
            config.totMaxNs = std::stod(argv[++index]);
        } else if (arg == "--tot-bins" && index + 1 < argc) {
            config.totBins = std::stoi(argv[++index]);
        } else if (arg == "--ch0-valid-rise-min-ns" && index + 1 < argc) {
            config.ch0ValidRiseMinNs = std::stod(argv[++index]);
        } else if (arg == "--ch0-valid-rise-max-ns" && index + 1 < argc) {
            config.ch0ValidRiseMaxNs = std::stod(argv[++index]);
        } else if (arg == "--ch0-valid-tot-min-ns" && index + 1 < argc) {
            config.ch0ValidTotMinNs = std::stod(argv[++index]);
        } else if (arg == "--ch0-valid-tot-max-ns" && index + 1 < argc) {
            config.ch0ValidTotMaxNs = std::stod(argv[++index]);
        } else if (arg == "--help") {
            std::cout
                << "Usage: dogma_cleaned_all_channel_time_tot --input <cleaned_pulses.tsv> --output-root <Results_root>\n"
                << "  [--run-key <name>] [--trigger-ref-dir-name Trigger_ref_ToT] [--ch0ref-dir-name Ch0_ref_TOT]\n"
                << "  [--trigger-ref-min-ns -6000] [--trigger-ref-max-ns 6000] [--trigger-ref-bins 2400]\n"
                << "  [--ch0ref-min-ns -6000] [--ch0ref-max-ns 6000] [--ch0ref-bins 2400]\n"
                << "  [--tot-min-ns 0] [--tot-max-ns 128] [--tot-bins 1280]\n"
                << "  [--ch0-valid-rise-min-ns -410] [--ch0-valid-rise-max-ns -395]\n"
                << "  [--ch0-valid-tot-min-ns 16.5] [--ch0-valid-tot-max-ns 19.5]\n";
            std::exit(0);
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
    if (!(config.triggerRefMinNs < config.triggerRefMaxNs)) {
        throw std::runtime_error("trigger-ref-min-ns must be smaller than trigger-ref-max-ns");
    }
    if (!(config.ch0RefMinNs < config.ch0RefMaxNs)) {
        throw std::runtime_error("ch0ref-min-ns must be smaller than ch0ref-max-ns");
    }
    if (!(config.totMinNs < config.totMaxNs)) {
        throw std::runtime_error("tot-min-ns must be smaller than tot-max-ns");
    }
    if (config.triggerRefBins <= 0 || config.ch0RefBins <= 0 || config.totBins <= 0) {
        throw std::runtime_error("Histogram bin counts must be positive");
    }
    if (!(config.ch0ValidRiseMinNs < config.ch0ValidRiseMaxNs)) {
        throw std::runtime_error("ch0-valid-rise-min-ns must be smaller than ch0-valid-rise-max-ns");
    }
    if (!(config.ch0ValidTotMinNs < config.ch0ValidTotMaxNs)) {
        throw std::runtime_error("ch0-valid-tot-min-ns must be smaller than ch0-valid-tot-max-ns");
    }
    return config;
}

bool parse_metadata_line(const std::string& line, std::string& key, std::string& value) {
    if (!starts_with(line, "#") || line.find('=') == std::string::npos) {
        return false;
    }
    const std::size_t separator = line.find('=');
    key = trim(std::string_view(line).substr(1, separator - 1));
    value = trim(std::string_view(line).substr(separator + 1));
    return true;
}

bool parse_pulse_row(const std::string& line, PulseRow& row) {
    std::istringstream stream(line);
    int tdcOrdinal = 0;
    int riseIndex = 0;
    int fallIndex = 0;
    double globalTriggerSeconds = 0.0;
    double fallNs = 0.0;
    if (!(stream >> row.windowIndex >> tdcOrdinal >> globalTriggerSeconds >> row.channel >> riseIndex >> fallIndex >> row.riseNs >> fallNs >> row.totNs)) {
        return false;
    }
    return row.channel >= 0 && row.channel < kChannelCount;
}

int histogram_bin_index(double value, int bins, double min, double max) {
    if (value < min) {
        return 0;
    }
    if (value >= max) {
        return bins + 1;
    }
    const double fraction = (value - min) / (max - min);
    return 1 + std::min(bins - 1, static_cast<int>(fraction * bins));
}

std::uint64_t sparse_key(int xBinZeroBased, int yBinZeroBased, int yBins) {
    return static_cast<std::uint64_t>(xBinZeroBased) * static_cast<std::uint64_t>(yBins)
         + static_cast<std::uint64_t>(yBinZeroBased);
}

void fill_sparse_histogram(SparseHistogram& histogram,
                           int channel,
                           int xBin,
                           int yBin,
                           int yBins) {
    if (xBin <= 0 || yBin <= 0) {
        return;
    }
    const std::uint64_t key = sparse_key(xBin - 1, yBin - 1, yBins);
    histogram[static_cast<std::size_t>(channel)][key] += 1;
}

void update_out_of_range_stats(std::uint64_t& outsideTime,
                               std::uint64_t& outsideTot,
                               std::uint64_t& outsideBoth,
                               int xBin,
                               int xBins,
                               int yBin,
                               int yBins) {
    const bool xOutside = xBin == 0 || xBin == xBins + 1;
    const bool yOutside = yBin == 0 || yBin == yBins + 1;
    if (xOutside && yOutside) {
        outsideBoth += 1;
    } else if (xOutside) {
        outsideTime += 1;
    } else if (yOutside) {
        outsideTot += 1;
    }
}

double tot_bin_width_ns(const Config& config) {
    return (config.totMaxNs - config.totMinNs) / static_cast<double>(config.totBins);
}

void fill_tot_distribution(const Config& config,
                           const PulseRow& row,
                           Stats& stats,
                           TotDistribution& totDistribution) {
    stats.observedTotMinNs = std::min(stats.observedTotMinNs, row.totNs);
    stats.observedTotMaxNs = std::max(stats.observedTotMaxNs, row.totNs);

    if (row.totNs < config.totMinNs) {
        stats.wholeRunTotRejectedBelowMin += 1;
        return;
    }

    const double binWidthNs = tot_bin_width_ns(config);
    const std::size_t totBin = static_cast<std::size_t>((row.totNs - config.totMinNs) / binWidthNs);
    auto& channelBins = totDistribution[static_cast<std::size_t>(row.channel)];
    if (totBin >= channelBins.size()) {
        channelBins.resize(totBin + 1, 0);
    }
    channelBins[totBin] += 1;
    stats.wholeRunTotAcceptedPulses += 1;
    stats.wholeRunTotAcceptedByChannel[static_cast<std::size_t>(row.channel)] += 1;
}

void fill_trigger_ref(const Config& config,
                      const PulseRow& row,
                      Stats& stats,
                      SparseHistogram& triggerRefHistogram) {
    const int xBin = histogram_bin_index(row.riseNs, config.triggerRefBins, config.triggerRefMinNs, config.triggerRefMaxNs);
    const int yBin = histogram_bin_index(row.totNs, config.totBins, config.totMinNs, config.totMaxNs);
    if (xBin > 0 && xBin <= config.triggerRefBins && yBin > 0 && yBin <= config.totBins) {
        fill_sparse_histogram(triggerRefHistogram, row.channel, xBin, yBin, config.totBins);
        stats.triggerRefAcceptedPulses += 1;
        stats.triggerRefAcceptedByChannel[static_cast<std::size_t>(row.channel)] += 1;
        return;
    }
    update_out_of_range_stats(stats.triggerRefRejectedOutsideTimeRange,
                              stats.triggerRefRejectedOutsideTotRange,
                              stats.triggerRefRejectedOutsideTimeAndTotRange,
                              xBin,
                              config.triggerRefBins,
                              yBin,
                              config.totBins);
}

void fill_ch0_ref_hit(const Config& config,
                      int channel,
                      double relativeTimeNs,
                      double totNs,
                      Stats& stats,
                      SparseHistogram& ch0RefHistogram) {
    const int xBin = histogram_bin_index(relativeTimeNs, config.ch0RefBins, config.ch0RefMinNs, config.ch0RefMaxNs);
    const int yBin = histogram_bin_index(totNs, config.totBins, config.totMinNs, config.totMaxNs);
    if (xBin > 0 && xBin <= config.ch0RefBins && yBin > 0 && yBin <= config.totBins) {
        fill_sparse_histogram(ch0RefHistogram, channel, xBin, yBin, config.totBins);
        stats.ch0RefAcceptedPulses += 1;
        stats.ch0RefAcceptedByChannel[static_cast<std::size_t>(channel)] += 1;
        return;
    }
    update_out_of_range_stats(stats.ch0RefRejectedOutsideTimeRange,
                              stats.ch0RefRejectedOutsideTotRange,
                              stats.ch0RefRejectedOutsideTimeAndTotRange,
                              xBin,
                              config.ch0RefBins,
                              yBin,
                              config.totBins);
}

void finalize_window(const WindowState& window,
                     const Config& config,
                     Stats& stats,
                     SparseHistogram& ch0RefHistogram) {
    if (!window.hasValue) {
        return;
    }

    stats.windowsSeen += 1;
    const auto& ch0Hits = window.hitsByChannel[0];
    if (ch0Hits.empty()) {
        stats.windowsWithoutCh0 += 1;
        return;
    }

    stats.totalCh0Candidates += ch0Hits.size();
    if (ch0Hits.size() > 1) {
        stats.windowsWithMultipleCh0 += 1;
    }

    const PulseHit* referenceHit = nullptr;
    std::size_t validCh0CandidateCount = 0;
    for (const PulseHit& hit : ch0Hits) {
        if (hit.riseNs < config.ch0ValidRiseMinNs || hit.riseNs > config.ch0ValidRiseMaxNs) {
            stats.rejectedCh0OutsideTimeRange += 1;
            continue;
        }
        if (hit.totNs < config.ch0ValidTotMinNs || hit.totNs > config.ch0ValidTotMaxNs) {
            stats.rejectedCh0OutsideTotRange += 1;
            continue;
        }
        stats.totalValidCh0Candidates += 1;
        validCh0CandidateCount += 1;
        if (referenceHit == nullptr) {
            referenceHit = &hit;
        }
    }

    if (referenceHit == nullptr) {
        stats.windowsWithoutValidCh0Candidate += 1;
        return;
    }

    stats.validCh0Windows += 1;
    if (validCh0CandidateCount == 1) {
        stats.windowsWithExactlyOneValidCh0Candidate += 1;
    } else {
        stats.windowsWithMultipleValidCh0Candidates += 1;
    }

    fill_ch0_ref_hit(config, 0, 0.0, referenceHit->totNs, stats, ch0RefHistogram);

    for (int channel = 1; channel < kChannelCount; ++channel) {
        for (const PulseHit& hit : window.hitsByChannel[static_cast<std::size_t>(channel)]) {
            fill_ch0_ref_hit(config, channel, hit.riseNs - referenceHit->riseNs, hit.totNs, stats, ch0RefHistogram);
        }
    }
}

void write_sparse_histogram(const fs::path& outputPath,
                            const std::string& referenceLabel,
                            double xMinNs,
                            double xMaxNs,
                            int xBins,
                            const Config& config,
                            const Stats& stats,
                            const SparseHistogram& histogram) {
    std::ofstream output(outputPath);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write sparse histogram: " + outputPath.string());
    }

    output << "# input_pulse_table=" << config.inputPath << '\n';
    if (!stats.sourceInputFile.empty()) {
        output << "# input_raw_file=" << stats.sourceInputFile << '\n';
    }
    output << "# run_key=" << config.runKey << '\n';
    output << "# x_reference=" << referenceLabel << '\n';
    output << "# x_min_ns=" << format_double(xMinNs) << '\n';
    output << "# x_max_ns=" << format_double(xMaxNs) << '\n';
    output << "# x_bins=" << xBins << '\n';
    output << "# tot_min_ns=" << format_double(config.totMinNs) << '\n';
    output << "# tot_max_ns=" << format_double(config.totMaxNs) << '\n';
    output << "# tot_bins=" << config.totBins << '\n';
    output << "# channel_count=" << kChannelCount << '\n';
    output << "# ch0_valid_rise_min_ns=" << format_double(config.ch0ValidRiseMinNs) << '\n';
    output << "# ch0_valid_rise_max_ns=" << format_double(config.ch0ValidRiseMaxNs) << '\n';
    output << "# ch0_valid_tot_min_ns=" << format_double(config.ch0ValidTotMinNs) << '\n';
    output << "# ch0_valid_tot_max_ns=" << format_double(config.ch0ValidTotMaxNs) << '\n';
    output << "# valid_ch0_windows=" << stats.validCh0Windows << '\n';
    output << "# windows_without_ch0=" << stats.windowsWithoutCh0 << '\n';
    output << "# windows_with_multiple_ch0=" << stats.windowsWithMultipleCh0 << '\n';
    output << "# windows_without_valid_ch0_candidate=" << stats.windowsWithoutValidCh0Candidate << '\n';
    output << "# windows_with_exactly_one_valid_ch0_candidate=" << stats.windowsWithExactlyOneValidCh0Candidate << '\n';
    output << "# windows_with_multiple_valid_ch0_candidates=" << stats.windowsWithMultipleValidCh0Candidates << '\n';
    output << "# total_ch0_candidates=" << stats.totalCh0Candidates << '\n';
    output << "# total_valid_ch0_candidates=" << stats.totalValidCh0Candidates << '\n';
    output << "# rejected_ch0_outside_time_range=" << stats.rejectedCh0OutsideTimeRange << '\n';
    output << "# rejected_ch0_outside_tot_range=" << stats.rejectedCh0OutsideTotRange << '\n';
    output << "# columns: channel x_bin_index_zero_based tot_bin_index_zero_based count\n";

    for (int channel = 0; channel < kChannelCount; ++channel) {
        std::vector<std::pair<std::uint64_t, std::uint64_t>> entries;
        entries.reserve(histogram[static_cast<std::size_t>(channel)].size());
        for (const auto& [key, count] : histogram[static_cast<std::size_t>(channel)]) {
            entries.emplace_back(key, count);
        }
        std::sort(entries.begin(), entries.end(), [](const auto& left, const auto& right) {
            return left.first < right.first;
        });
        for (const auto& [key, count] : entries) {
            const std::uint64_t xBin = key / static_cast<std::uint64_t>(config.totBins);
            const std::uint64_t yBin = key % static_cast<std::uint64_t>(config.totBins);
            output << channel << '\t' << xBin << '\t' << yBin << '\t' << count << '\n';
        }
    }
}

void write_tot_distribution(const fs::path& outputPath,
                            const Config& config,
                            const Stats& stats,
                            const TotDistribution& distribution) {
    std::ofstream output(outputPath);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write TOT distribution: " + outputPath.string());
    }

    const double binWidthNs = tot_bin_width_ns(config);
    std::size_t dynamicBins = 0;
    for (int channel = 0; channel < kChannelCount; ++channel) {
        dynamicBins = std::max(dynamicBins, distribution[static_cast<std::size_t>(channel)].size());
    }
    if (dynamicBins == 0) {
        dynamicBins = 1;
    }
    const double dynamicTotMaxNs = config.totMinNs + binWidthNs * static_cast<double>(dynamicBins);

    output << "# input_pulse_table=" << config.inputPath << '\n';
    if (!stats.sourceInputFile.empty()) {
        output << "# input_raw_file=" << stats.sourceInputFile << '\n';
    }
    output << "# run_key=" << config.runKey << '\n';
    output << "# distribution_source=whole_run_parsed_pulses" << '\n';
    output << "# tot_min_ns=" << format_double(config.totMinNs) << '\n';
    output << "# tot_max_ns=" << format_double(dynamicTotMaxNs) << '\n';
    output << "# tot_bins=" << dynamicBins << '\n';
    output << "# tot_bin_width_ns=" << format_double(binWidthNs) << '\n';
    output << "# channel_count=" << kChannelCount << '\n';
    output << "# whole_run_tot_accepted_pulses=" << stats.wholeRunTotAcceptedPulses << '\n';
    output << "# whole_run_tot_rejected_below_min=" << stats.wholeRunTotRejectedBelowMin << '\n';
    if (stats.parsedPulseRows > 0) {
        output << "# observed_tot_min_ns=" << format_double(stats.observedTotMinNs) << '\n';
        output << "# observed_tot_max_ns=" << format_double(stats.observedTotMaxNs) << '\n';
    }
    output << "# columns: channel tot_bin_index_zero_based count\n";

    for (int channel = 0; channel < kChannelCount; ++channel) {
        const auto& channelBins = distribution[static_cast<std::size_t>(channel)];
        for (std::size_t totBin = 0; totBin < channelBins.size(); ++totBin) {
            const std::uint64_t count = channelBins[totBin];
            if (count == 0) {
                continue;
            }
            output << channel << '\t' << totBin << '\t' << count << '\n';
        }
    }
}

void write_summary(const fs::path& outputPath,
                   const Config& config,
                   const Stats& stats,
                   const fs::path& triggerRefPath,
                   const fs::path& ch0RefPath,
                   const fs::path& totDistribPath) {
    std::ofstream output(outputPath);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write cleaned TOT summary: " + outputPath.string());
    }

    output << "input_pulse_table=" << config.inputPath << '\n';
    if (!stats.sourceInputFile.empty()) {
        output << "input_raw_file=" << stats.sourceInputFile << '\n';
    }
    output << "run_key=" << config.runKey << '\n';
    output << "trigger_ref_output=" << triggerRefPath << '\n';
    output << "ch0_ref_output=" << ch0RefPath << '\n';
    output << "tot_distrib_output=" << totDistribPath << '\n';
    output << "trigger_ref_min_ns=" << format_double(config.triggerRefMinNs) << '\n';
    output << "trigger_ref_max_ns=" << format_double(config.triggerRefMaxNs) << '\n';
    output << "trigger_ref_bins=" << config.triggerRefBins << '\n';
    output << "ch0_ref_min_ns=" << format_double(config.ch0RefMinNs) << '\n';
    output << "ch0_ref_max_ns=" << format_double(config.ch0RefMaxNs) << '\n';
    output << "ch0_ref_bins=" << config.ch0RefBins << '\n';
    output << "tot_min_ns=" << format_double(config.totMinNs) << '\n';
    output << "tot_max_ns=" << format_double(config.totMaxNs) << '\n';
    output << "tot_bins=" << config.totBins << '\n';
    output << "ch0_valid_rise_min_ns=" << format_double(config.ch0ValidRiseMinNs) << '\n';
    output << "ch0_valid_rise_max_ns=" << format_double(config.ch0ValidRiseMaxNs) << '\n';
    output << "ch0_valid_tot_min_ns=" << format_double(config.ch0ValidTotMinNs) << '\n';
    output << "ch0_valid_tot_max_ns=" << format_double(config.ch0ValidTotMaxNs) << '\n';
    if (stats.parsedPulseRows > 0) {
        output << "observed_tot_min_ns=" << format_double(stats.observedTotMinNs) << '\n';
        output << "observed_tot_max_ns=" << format_double(stats.observedTotMaxNs) << '\n';
    }
    output << "parsed_pulse_rows=" << stats.parsedPulseRows << '\n';
    output << "invalid_rows=" << stats.invalidRows << '\n';
    output << "windows_seen=" << stats.windowsSeen << '\n';
    output << "valid_ch0_windows=" << stats.validCh0Windows << '\n';
    output << "windows_without_ch0=" << stats.windowsWithoutCh0 << '\n';
    output << "windows_with_multiple_ch0=" << stats.windowsWithMultipleCh0 << '\n';
    output << "windows_without_valid_ch0_candidate=" << stats.windowsWithoutValidCh0Candidate << '\n';
    output << "windows_with_exactly_one_valid_ch0_candidate=" << stats.windowsWithExactlyOneValidCh0Candidate << '\n';
    output << "windows_with_multiple_valid_ch0_candidates=" << stats.windowsWithMultipleValidCh0Candidates << '\n';
    output << "total_ch0_candidates=" << stats.totalCh0Candidates << '\n';
    output << "total_valid_ch0_candidates=" << stats.totalValidCh0Candidates << '\n';
    output << "rejected_ch0_outside_time_range=" << stats.rejectedCh0OutsideTimeRange << '\n';
    output << "rejected_ch0_outside_tot_range=" << stats.rejectedCh0OutsideTotRange << '\n';
    output << "trigger_ref_accepted_pulses=" << stats.triggerRefAcceptedPulses << '\n';
    output << "trigger_ref_rejected_outside_time_range=" << stats.triggerRefRejectedOutsideTimeRange << '\n';
    output << "trigger_ref_rejected_outside_tot_range=" << stats.triggerRefRejectedOutsideTotRange << '\n';
    output << "trigger_ref_rejected_outside_time_and_tot_range=" << stats.triggerRefRejectedOutsideTimeAndTotRange << '\n';
    output << "ch0_ref_accepted_pulses=" << stats.ch0RefAcceptedPulses << '\n';
    output << "ch0_ref_rejected_outside_time_range=" << stats.ch0RefRejectedOutsideTimeRange << '\n';
    output << "ch0_ref_rejected_outside_tot_range=" << stats.ch0RefRejectedOutsideTotRange << '\n';
    output << "ch0_ref_rejected_outside_time_and_tot_range=" << stats.ch0RefRejectedOutsideTimeAndTotRange << '\n';
        output << "whole_run_tot_accepted_pulses=" << stats.wholeRunTotAcceptedPulses << '\n';
        output << "whole_run_tot_rejected_below_min=" << stats.wholeRunTotRejectedBelowMin << '\n';

    for (int channel = 0; channel < kChannelCount; ++channel) {
        output << "channel_" << std::setfill('0') << std::setw(2) << channel << "_pulses="
               << stats.pulseCountByChannel[static_cast<std::size_t>(channel)] << '\n';
        output << "channel_" << std::setfill('0') << std::setw(2) << channel << "_trigger_ref_accepted="
               << stats.triggerRefAcceptedByChannel[static_cast<std::size_t>(channel)] << '\n';
        output << "channel_" << std::setfill('0') << std::setw(2) << channel << "_ch0_ref_accepted="
               << stats.ch0RefAcceptedByChannel[static_cast<std::size_t>(channel)] << '\n';
         output << "channel_" << std::setfill('0') << std::setw(2) << channel << "_whole_run_tot_accepted="
             << stats.wholeRunTotAcceptedByChannel[static_cast<std::size_t>(channel)] << '\n';
    }
    output << std::setfill(' ');
}

}  // namespace

int main(int argc, char** argv) {
    try {
        Config config = parse_args(argc, argv);
        std::ifstream input(config.inputPath);
        if (!input.is_open()) {
            throw std::runtime_error("Unable to open cleaned pulse table: " + config.inputPath);
        }

        Stats stats;
        SparseHistogram triggerRefHistogram;
        SparseHistogram ch0RefHistogram;
        TotDistribution totDistribution;
        WindowState currentWindow;

        std::string line;
        while (std::getline(input, line)) {
            if (line.empty()) {
                continue;
            }

            std::string key;
            std::string value;
            if (parse_metadata_line(line, key, value)) {
                if (key == "run_key" && config.runKey.empty()) {
                    config.runKey = value;
                } else if (key == "input_file") {
                    stats.sourceInputFile = value;
                }
                continue;
            }
            if (starts_with(line, "window_index")) {
                continue;
            }

            PulseRow row;
            if (!parse_pulse_row(line, row)) {
                stats.invalidRows += 1;
                continue;
            }

            stats.parsedPulseRows += 1;
            stats.pulseCountByChannel[static_cast<std::size_t>(row.channel)] += 1;
            fill_tot_distribution(config, row, stats, totDistribution);
            fill_trigger_ref(config, row, stats, triggerRefHistogram);

            if (!currentWindow.hasValue) {
                currentWindow.hasValue = true;
                currentWindow.windowIndex = row.windowIndex;
            }
            if (row.windowIndex != currentWindow.windowIndex) {
                finalize_window(currentWindow, config, stats, ch0RefHistogram);
                currentWindow = WindowState{};
                currentWindow.hasValue = true;
                currentWindow.windowIndex = row.windowIndex;
            }
            currentWindow.hitsByChannel[static_cast<std::size_t>(row.channel)].push_back(PulseHit{row.riseNs, row.totNs});
        }

        finalize_window(currentWindow, config, stats, ch0RefHistogram);

        if (config.runKey.empty()) {
            config.runKey = strip_known_suffixes(fs::path(config.inputPath).filename().string());
        }

        const fs::path runDir = fs::path(config.outputRoot) / config.runKey;
        const fs::path triggerRefDir = runDir / config.triggerRefDirName;
        const fs::path ch0RefDir = runDir / config.ch0RefDirName;
        const fs::path totDistribDir = runDir / config.totDistribDirName;
        fs::create_directories(triggerRefDir);
        fs::create_directories(ch0RefDir);
        fs::create_directories(totDistribDir);

        const fs::path triggerRefPath = triggerRefDir / (config.runKey + "_Trigger_ref_ToT_sparse.tsv");
        const fs::path ch0RefPath = ch0RefDir / (config.runKey + "_Ch0_ref_TOT_sparse.tsv");
        const fs::path totDistribPath = totDistribDir / (config.runKey + "_TOT_distrib.tsv");
        const fs::path summaryPath = runDir / (config.runKey + "_cleaned_time_tot_summary.txt");

        write_sparse_histogram(triggerRefPath,
                               "trigger_window_rise_time_ns",
                               config.triggerRefMinNs,
                               config.triggerRefMaxNs,
                               config.triggerRefBins,
                               config,
                               stats,
                               triggerRefHistogram);
        write_sparse_histogram(ch0RefPath,
                               "validated_first_ch0_reference_hit",
                               config.ch0RefMinNs,
                               config.ch0RefMaxNs,
                               config.ch0RefBins,
                               config,
                               stats,
                               ch0RefHistogram);
        write_tot_distribution(totDistribPath, config, stats, totDistribution);
        write_summary(summaryPath, config, stats, triggerRefPath, ch0RefPath, totDistribPath);

        std::cout << "Trigger-ref sparse histogram: " << triggerRefPath << '\n';
        std::cout << "Ch0-ref sparse histogram: " << ch0RefPath << '\n';
        std::cout << "ToT distribution histogram: " << totDistribPath << '\n';
        std::cout << "Summary: " << summaryPath << '\n';
        std::cout << "Parsed pulse rows: " << stats.parsedPulseRows << '\n';
        std::cout << "Valid ch0 windows: " << stats.validCh0Windows << '\n';
        std::cout << "Windows without valid ch0 candidate: " << stats.windowsWithoutValidCh0Candidate << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}