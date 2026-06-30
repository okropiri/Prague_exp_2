#include <algorithm>
#include <cmath>
#include <array>
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
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr int kChannelCount = 32;

struct Config {
    std::string inputPath;
    std::string outputRoot;
    std::string runKey;
    std::string absDirName = "Abs_rates";
    std::string ch0RefDirName = "Ch0_ref_Rates";
    double absDisplayBinWidthSeconds = 0.1;
    double ch0RefMinNs = -6000.0;
    double ch0RefMaxNs = 6000.0;
    int ch0RefBins = 12000;
    double ch0ValidRiseMinNs = -410.0;
    double ch0ValidRiseMaxNs = -395.0;
    double ch0ValidTotMinNs = 16.5;
    double ch0ValidTotMaxNs = 19.5;
};

struct PulseRow {
    std::uint64_t windowIndex = 0;
    double globalTriggerSeconds = 0.0;
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
    std::uint64_t negativeAbsoluteTimeRows = 0;
    bool hasFirstAbsoluteNs = false;
    std::int64_t firstAbsoluteNs = 0;
    std::int64_t lastAbsoluteNs = 0;
    std::array<std::uint64_t, kChannelCount> pulseCountByChannel{};
};

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
        } else if (arg == "--abs-dir-name" && index + 1 < argc) {
            config.absDirName = argv[++index];
        } else if (arg == "--ch0ref-dir-name" && index + 1 < argc) {
            config.ch0RefDirName = argv[++index];
        } else if (arg == "--abs-display-bin-width-seconds" && index + 1 < argc) {
            config.absDisplayBinWidthSeconds = std::stod(argv[++index]);
        } else if (arg == "--ch0ref-min-ns" && index + 1 < argc) {
            config.ch0RefMinNs = std::stod(argv[++index]);
        } else if (arg == "--ch0ref-max-ns" && index + 1 < argc) {
            config.ch0RefMaxNs = std::stod(argv[++index]);
        } else if (arg == "--ch0ref-bins" && index + 1 < argc) {
            config.ch0RefBins = std::stoi(argv[++index]);
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
                << "Usage: dogma_cleaned_all_channel_rates --input <cleaned_pulses.tsv> --output-root <Results_root>\n"
                << "  [--run-key <name>] [--abs-dir-name Abs_rates] [--ch0ref-dir-name Ch0_ref_Rates]\n"
                << "  [--abs-display-bin-width-seconds 0.1]\n"
                << "  [--ch0ref-min-ns -6000] [--ch0ref-max-ns 6000] [--ch0ref-bins 12000]\n"
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
    if (config.absDisplayBinWidthSeconds <= 0.0) {
        throw std::runtime_error("abs-display-bin-width-seconds must be positive");
    }
    if (!(config.ch0RefMinNs < config.ch0RefMaxNs)) {
        throw std::runtime_error("ch0ref-min-ns must be smaller than ch0ref-max-ns");
    }
    if (config.ch0RefBins <= 0) {
        throw std::runtime_error("ch0ref-bins must be positive");
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
    double fallNs = 0.0;
    if (!(stream >> row.windowIndex >> tdcOrdinal >> row.globalTriggerSeconds >> row.channel >> riseIndex >> fallIndex >> row.riseNs >> fallNs >> row.totNs)) {
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
    return 1 + std::min(bins - 1, static_cast<int>(std::floor(fraction * bins)));
}

void ensure_abs_capacity(std::array<std::vector<std::uint64_t>, kChannelCount>& byChannel,
                         std::vector<std::uint64_t>& total,
                         std::size_t size) {
    for (auto& values : byChannel) {
        if (values.size() < size) {
            values.resize(size, 0);
        }
    }
    if (total.size() < size) {
        total.resize(size, 0);
    }
}

void finalize_window(const WindowState& window,
                     const Config& config,
                     Stats& stats,
                     std::array<std::vector<std::uint64_t>, kChannelCount>& ch0RefCounts,
                     std::vector<std::uint64_t>& ch0RefTotals) {
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
    const double referenceNs = referenceHit->riseNs;

    for (int channel = 0; channel < kChannelCount; ++channel) {
        if (channel == 0) {
            const int binIndex = histogram_bin_index(0.0, config.ch0RefBins, config.ch0RefMinNs, config.ch0RefMaxNs);
            ch0RefCounts[0][static_cast<std::size_t>(binIndex)] += 1;
            ch0RefTotals[static_cast<std::size_t>(binIndex)] += 1;
            continue;
        }

        for (const PulseHit& hit : window.hitsByChannel[static_cast<std::size_t>(channel)]) {
            const int binIndex = histogram_bin_index(hit.riseNs - referenceNs, config.ch0RefBins, config.ch0RefMinNs, config.ch0RefMaxNs);
            ch0RefCounts[static_cast<std::size_t>(channel)][static_cast<std::size_t>(binIndex)] += 1;
            ch0RefTotals[static_cast<std::size_t>(binIndex)] += 1;
        }
    }
}

void write_abs_matrix(const fs::path& outputPath,
                      const Config& config,
                      const Stats& stats,
                      const std::array<std::vector<std::uint64_t>, kChannelCount>& absCounts,
                      const std::vector<std::uint64_t>& absTotals) {
    std::ofstream output(outputPath);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write absolute-rate matrix: " + outputPath.string());
    }

    output << "# input_pulse_table=" << config.inputPath << '\n';
    if (!stats.sourceInputFile.empty()) {
        output << "# input_raw_file=" << stats.sourceInputFile << '\n';
    }
    output << "# run_key=" << config.runKey << '\n';
    output << "# display_bin_width_seconds=" << format_double(config.absDisplayBinWidthSeconds) << '\n';
    output << "# run_start_absolute_seconds=" << format_double(static_cast<double>(stats.firstAbsoluteNs) * 1e-9) << '\n';
    output << "# last_absolute_seconds=" << format_double(static_cast<double>(stats.lastAbsoluteNs) * 1e-9) << '\n';
    output << "# parsed_pulse_rows=" << stats.parsedPulseRows << '\n';
    output << "# channel_count=" << kChannelCount << '\n';
    output << "# columns: bin_index time_center_seconds";
    for (int channel = 0; channel < kChannelCount; ++channel) {
        output << " ch" << std::setfill('0') << std::setw(2) << channel;
    }
    output << " total\n";
    output << std::setfill(' ');

    for (std::size_t index = 0; index < absTotals.size(); ++index) {
        const double timeCenterSeconds = (static_cast<double>(index) + 0.5) * config.absDisplayBinWidthSeconds;
        output << index << '\t' << format_double(timeCenterSeconds);
        for (int channel = 0; channel < kChannelCount; ++channel) {
            output << '\t' << absCounts[static_cast<std::size_t>(channel)][index];
        }
        output << '\t' << absTotals[index] << '\n';
    }
}

void write_ch0ref_matrix(const fs::path& outputPath,
                         const Config& config,
                         const Stats& stats,
                         const std::array<std::vector<std::uint64_t>, kChannelCount>& ch0RefCounts,
                         const std::vector<std::uint64_t>& ch0RefTotals) {
    std::ofstream output(outputPath);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write ch0-referenced matrix: " + outputPath.string());
    }

    output << "# input_pulse_table=" << config.inputPath << '\n';
    if (!stats.sourceInputFile.empty()) {
        output << "# input_raw_file=" << stats.sourceInputFile << '\n';
    }
    output << "# run_key=" << config.runKey << '\n';
    output << "# x_min_ns=" << format_double(config.ch0RefMinNs) << '\n';
    output << "# x_max_ns=" << format_double(config.ch0RefMaxNs) << '\n';
    output << "# bins=" << config.ch0RefBins << '\n';
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
    output << "# channel_count=" << kChannelCount << '\n';
    output << "# columns: bin_index time_center_ns";
    for (int channel = 0; channel < kChannelCount; ++channel) {
        output << " ch" << std::setfill('0') << std::setw(2) << channel;
    }
    output << " total\n";
    output << std::setfill(' ');

    const double binWidthNs = (config.ch0RefMaxNs - config.ch0RefMinNs) / static_cast<double>(config.ch0RefBins);
    for (int bin = 1; bin <= config.ch0RefBins; ++bin) {
        const double centerNs = config.ch0RefMinNs + (static_cast<double>(bin) - 0.5) * binWidthNs;
        output << bin << '\t' << format_double(centerNs);
        for (int channel = 0; channel < kChannelCount; ++channel) {
            output << '\t' << ch0RefCounts[static_cast<std::size_t>(channel)][static_cast<std::size_t>(bin)];
        }
        output << '\t' << ch0RefTotals[static_cast<std::size_t>(bin)] << '\n';
    }
}

void write_summary(const fs::path& outputPath, const Config& config, const Stats& stats) {
    std::ofstream output(outputPath);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write cleaned-rate summary: " + outputPath.string());
    }

    output << "input_pulse_table=" << config.inputPath << '\n';
    if (!stats.sourceInputFile.empty()) {
        output << "input_raw_file=" << stats.sourceInputFile << '\n';
    }
    output << "run_key=" << config.runKey << '\n';
    output << "abs_display_bin_width_seconds=" << format_double(config.absDisplayBinWidthSeconds) << '\n';
    output << "ch0ref_min_ns=" << format_double(config.ch0RefMinNs) << '\n';
    output << "ch0ref_max_ns=" << format_double(config.ch0RefMaxNs) << '\n';
    output << "ch0ref_bins=" << config.ch0RefBins << '\n';
    output << "ch0_valid_rise_min_ns=" << format_double(config.ch0ValidRiseMinNs) << '\n';
    output << "ch0_valid_rise_max_ns=" << format_double(config.ch0ValidRiseMaxNs) << '\n';
    output << "ch0_valid_tot_min_ns=" << format_double(config.ch0ValidTotMinNs) << '\n';
    output << "ch0_valid_tot_max_ns=" << format_double(config.ch0ValidTotMaxNs) << '\n';
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
    output << "negative_absolute_time_rows=" << stats.negativeAbsoluteTimeRows << '\n';
    if (stats.hasFirstAbsoluteNs) {
        output << "run_start_absolute_seconds=" << format_double(static_cast<double>(stats.firstAbsoluteNs) * 1e-9) << '\n';
        output << "last_absolute_seconds=" << format_double(static_cast<double>(stats.lastAbsoluteNs) * 1e-9) << '\n';
        output << "covered_duration_seconds=" << format_double(static_cast<double>(stats.lastAbsoluteNs - stats.firstAbsoluteNs) * 1e-9) << '\n';
    }
    for (int channel = 0; channel < kChannelCount; ++channel) {
        output << "channel_" << std::setfill('0') << std::setw(2) << channel << "_pulses="
               << stats.pulseCountByChannel[static_cast<std::size_t>(channel)] << '\n';
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
        std::array<std::vector<std::uint64_t>, kChannelCount> absCounts;
        std::vector<std::uint64_t> absTotals;
        std::array<std::vector<std::uint64_t>, kChannelCount> ch0RefCounts;
        std::vector<std::uint64_t> ch0RefTotals(static_cast<std::size_t>(config.ch0RefBins + 2), 0);
        for (auto& values : ch0RefCounts) {
            values.assign(static_cast<std::size_t>(config.ch0RefBins + 2), 0);
        }

        WindowState currentWindow;
        const std::int64_t absDisplayBinWidthNs = static_cast<std::int64_t>(std::llround(config.absDisplayBinWidthSeconds * 1e9));
        if (absDisplayBinWidthNs <= 0) {
            throw std::runtime_error("abs-display-bin-width-seconds is too small");
        }

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

            const std::int64_t absoluteNs = static_cast<std::int64_t>(std::llround(row.globalTriggerSeconds * 1e9 + row.riseNs));
            if (!stats.hasFirstAbsoluteNs) {
                stats.hasFirstAbsoluteNs = true;
                stats.firstAbsoluteNs = absoluteNs;
                stats.lastAbsoluteNs = absoluteNs;
            }
            if (absoluteNs < stats.firstAbsoluteNs) {
                stats.negativeAbsoluteTimeRows += 1;
                continue;
            }
            stats.lastAbsoluteNs = std::max(stats.lastAbsoluteNs, absoluteNs);
            const std::size_t absIndex = static_cast<std::size_t>((absoluteNs - stats.firstAbsoluteNs) / absDisplayBinWidthNs);
            ensure_abs_capacity(absCounts, absTotals, absIndex + 1);
            absCounts[static_cast<std::size_t>(row.channel)][absIndex] += 1;
            absTotals[absIndex] += 1;

            if (!currentWindow.hasValue) {
                currentWindow.hasValue = true;
                currentWindow.windowIndex = row.windowIndex;
            }
            if (row.windowIndex != currentWindow.windowIndex) {
                finalize_window(currentWindow, config, stats, ch0RefCounts, ch0RefTotals);
                currentWindow = WindowState{};
                currentWindow.hasValue = true;
                currentWindow.windowIndex = row.windowIndex;
            }
            currentWindow.hitsByChannel[static_cast<std::size_t>(row.channel)].push_back(PulseHit{row.riseNs, row.totNs});
        }

        finalize_window(currentWindow, config, stats, ch0RefCounts, ch0RefTotals);

        if (config.runKey.empty()) {
            config.runKey = strip_known_suffixes(fs::path(config.inputPath).filename().string());
        }

        const fs::path runDir = fs::path(config.outputRoot) / config.runKey;
        const fs::path absDir = runDir / config.absDirName;
        const fs::path ch0RefDir = runDir / config.ch0RefDirName;
        fs::create_directories(absDir);
        fs::create_directories(ch0RefDir);

        const fs::path absMatrixPath = absDir / (config.runKey + "_Abs_rates_matrix.tsv");
        const fs::path ch0RefMatrixPath = ch0RefDir / (config.runKey + "_Ch0_ref_Rates_matrix.tsv");
        const fs::path summaryPath = runDir / (config.runKey + "_cleaned_rates_summary.txt");

        write_abs_matrix(absMatrixPath, config, stats, absCounts, absTotals);
        write_ch0ref_matrix(ch0RefMatrixPath, config, stats, ch0RefCounts, ch0RefTotals);
        write_summary(summaryPath, config, stats);

        std::cout << "Absolute matrix: " << absMatrixPath << '\n';
        std::cout << "Ch0-ref matrix: " << ch0RefMatrixPath << '\n';
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