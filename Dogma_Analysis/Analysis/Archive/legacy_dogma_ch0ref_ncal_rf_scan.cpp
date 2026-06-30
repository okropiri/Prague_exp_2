#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <deque>
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
#include <utility>
#include <vector>

namespace {

struct Config {
    std::string inputPath;
    std::string outputPrefix;
    double windowMinNs = -6000.0;
    double windowMaxNs = 6000.0;
    double totMinNs = 0.0;
    double totMaxNs = 128.0;
    int totBins = 512;
    double scoreTotMinNs = 0.0;
    double scoreTotMaxNs = 128.0;
    double initialPeriodMinNs = 39.0;
    double initialPeriodMaxNs = 41.0;
    double initialStepNs = 0.05;
    int refineRounds = 3;
    int refineHalfSpanSteps = 5;
    double refineFactor = 10.0;
    double phaseBinWidthNs = 0.25;
    double peakWindowNs = 8.0;
    std::uint64_t minSelectedPulses = 1000;
    std::uint64_t scoreStride = 1;
    double triggerResetThresholdSeconds = 1.0;
};

struct PulsePair {
    double riseNs = 0.0;
    double totNs = 0.0;
};

struct StoredPulse {
    float timeNs = 0.0F;
    float totNs = 0.0F;
};

struct ScoredPulse {
    float timeNs = 0.0F;
    double globalTimeSeconds = 0.0;
};

struct WindowState {
    std::vector<double> ch0Rises;
    std::deque<double> pendingNcal1Rises;
    std::vector<PulsePair> ncal1Pairs;
    bool hasGlobalTriggerSeconds = false;
    double globalTriggerSeconds = 0.0;
};

struct AnalysisStats {
    std::uint64_t parsedHeaders = 0;
    std::uint64_t parsedRows = 0;
    std::uint64_t finalizedWindows = 0;
    std::uint64_t validCh0Rises = 0;
    std::uint64_t validNcal1Rises = 0;
    std::uint64_t validNcal1Pairs = 0;
    std::uint64_t windowsWithValidCh0 = 0;
    std::uint64_t windowsWithoutCh0 = 0;
    std::uint64_t windowsWithMultipleCh0 = 0;
    bool duplicateRestartDetected = false;
    double duplicateRestartPreviousTriggerSeconds = 0.0;
    double duplicateRestartNewTriggerSeconds = 0.0;
    std::uint64_t duplicateRestartHeadersSkipped = 0;
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

struct CandidateMetrics {
    int roundIndex = 0;
    double periodNs = 0.0;
    double phaseOriginNs = 0.0;
    double peakCenterNs = 0.0;
    std::uint64_t peakHeight = 0;
    std::uint64_t selectedPulses = 0;
    double selectedFraction = 0.0;
    double sigmaNs = std::numeric_limits<double>::infinity();
    double meanResidualNs = 0.0;
    double coherence = 0.0;
    double driftSlopeNsPerCycle = 0.0;
    double driftInterceptNs = 0.0;
    double merit = 0.0;
    bool valid = false;
};

struct CycleResidualStats {
    std::int64_t cycleIndex = 0;
    std::uint64_t count = 0;
    double meanResidualNs = 0.0;
    double rmsResidualNs = 0.0;
};

struct CycleResidualPoint {
    std::int64_t cycleIndex = 0;
    double residualNs = 0.0;
    double globalTimeSeconds = 0.0;
};

constexpr int kTriggerOrdinal = 0;
constexpr int kSignalOrdinal = 1;
constexpr int kCh0Channel = 0;
constexpr int kNcal1Channel = 2;
constexpr double kMaxAcceptedHitTimeNs = 50000.0;
constexpr double kTwoPi = 6.28318530717958647692;

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
        } else if (arg == "--window-min-ns" && index + 1 < argc) {
            config.windowMinNs = std::stod(argv[++index]);
        } else if (arg == "--window-max-ns" && index + 1 < argc) {
            config.windowMaxNs = std::stod(argv[++index]);
        } else if (arg == "--tot-min-ns" && index + 1 < argc) {
            config.totMinNs = std::stod(argv[++index]);
        } else if (arg == "--tot-max-ns" && index + 1 < argc) {
            config.totMaxNs = std::stod(argv[++index]);
        } else if (arg == "--tot-bins" && index + 1 < argc) {
            config.totBins = std::stoi(argv[++index]);
        } else if (arg == "--score-tot-min-ns" && index + 1 < argc) {
            config.scoreTotMinNs = std::stod(argv[++index]);
        } else if (arg == "--score-tot-max-ns" && index + 1 < argc) {
            config.scoreTotMaxNs = std::stod(argv[++index]);
        } else if (arg == "--initial-period-min-ns" && index + 1 < argc) {
            config.initialPeriodMinNs = std::stod(argv[++index]);
        } else if (arg == "--initial-period-max-ns" && index + 1 < argc) {
            config.initialPeriodMaxNs = std::stod(argv[++index]);
        } else if (arg == "--initial-step-ns" && index + 1 < argc) {
            config.initialStepNs = std::stod(argv[++index]);
        } else if (arg == "--refine-rounds" && index + 1 < argc) {
            config.refineRounds = std::stoi(argv[++index]);
        } else if (arg == "--refine-half-span-steps" && index + 1 < argc) {
            config.refineHalfSpanSteps = std::stoi(argv[++index]);
        } else if (arg == "--refine-factor" && index + 1 < argc) {
            config.refineFactor = std::stod(argv[++index]);
        } else if (arg == "--phase-bin-width-ns" && index + 1 < argc) {
            config.phaseBinWidthNs = std::stod(argv[++index]);
        } else if (arg == "--peak-window-ns" && index + 1 < argc) {
            config.peakWindowNs = std::stod(argv[++index]);
        } else if (arg == "--min-selected-pulses" && index + 1 < argc) {
            config.minSelectedPulses = static_cast<std::uint64_t>(std::stoull(argv[++index]));
        } else if (arg == "--score-stride" && index + 1 < argc) {
            config.scoreStride = static_cast<std::uint64_t>(std::stoull(argv[++index]));
        } else if (arg == "--trigger-reset-threshold-seconds" && index + 1 < argc) {
            config.triggerResetThresholdSeconds = std::stod(argv[++index]);
        } else if (arg == "--help") {
            std::cout
                << "Usage: dogma_ch0ref_ncal_rf_scan --input <file> --output-prefix <prefix>\n"
                << "  [--window-min-ns -6000] [--window-max-ns 6000]\n"
                << "  [--tot-min-ns 0] [--tot-max-ns 128] [--tot-bins 512]\n"
                << "  [--score-tot-min-ns 0] [--score-tot-max-ns 128]\n"
                << "  [--initial-period-min-ns 39.0] [--initial-period-max-ns 41.0] [--initial-step-ns 0.05]\n"
                << "  [--refine-rounds 3] [--refine-half-span-steps 5] [--refine-factor 10]\n"
                << "  [--phase-bin-width-ns 0.25] [--peak-window-ns 8.0] [--min-selected-pulses 1000] [--score-stride 1]\n"
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
    if (!(config.windowMinNs < config.windowMaxNs)) {
        throw std::runtime_error("window-min-ns must be smaller than window-max-ns");
    }
    if (!(config.totMinNs < config.totMaxNs)) {
        throw std::runtime_error("tot-min-ns must be smaller than tot-max-ns");
    }
    if (config.totBins <= 0) {
        throw std::runtime_error("tot-bins must be positive");
    }
    if (!(config.scoreTotMinNs < config.scoreTotMaxNs)) {
        throw std::runtime_error("score-tot-min-ns must be smaller than score-tot-max-ns");
    }
    if (!(config.initialPeriodMinNs < config.initialPeriodMaxNs)) {
        throw std::runtime_error("initial-period-min-ns must be smaller than initial-period-max-ns");
    }
    if (config.initialStepNs <= 0.0) {
        throw std::runtime_error("initial-step-ns must be positive");
    }
    if (config.refineRounds <= 0) {
        throw std::runtime_error("refine-rounds must be positive");
    }
    if (config.refineHalfSpanSteps <= 0) {
        throw std::runtime_error("refine-half-span-steps must be positive");
    }
    if (config.refineFactor <= 1.0) {
        throw std::runtime_error("refine-factor must be larger than 1");
    }
    if (config.phaseBinWidthNs <= 0.0) {
        throw std::runtime_error("phase-bin-width-ns must be positive");
    }
    if (config.peakWindowNs <= 0.0) {
        throw std::runtime_error("peak-window-ns must be positive");
    }
    if (config.triggerResetThresholdSeconds <= 0.0) {
        throw std::runtime_error("trigger-reset-threshold-seconds must be positive");
    }
    if (config.scoreStride == 0) {
        throw std::runtime_error("score-stride must be positive");
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

double positive_mod(double value, double period) {
    const double wrapped = std::fmod(value, period);
    return wrapped < 0.0 ? wrapped + period : wrapped;
}

double wrap_centered(double value, double period) {
    return positive_mod(value + 0.5 * period, period) - 0.5 * period;
}

bool passes_score_tot_gate(const Config& config, float totNs) {
    return static_cast<double>(totNs) >= config.scoreTotMinNs && static_cast<double>(totNs) < config.scoreTotMaxNs;
}

void finalize_window(const WindowState& window,
                     const Config& config,
                     AnalysisStats& stats,
                     std::vector<StoredPulse>& pulses,
                     std::vector<ScoredPulse>& scorePulses,
                     std::uint64_t& scoreOrdinal) {
    stats.finalizedWindows += 1;
    if (window.ch0Rises.empty()) {
        stats.windowsWithoutCh0 += 1;
        return;
    }

    stats.windowsWithValidCh0 += 1;
    if (window.ch0Rises.size() > 1) {
        stats.windowsWithMultipleCh0 += 1;
    }

    const double referenceTimeNs = *std::min_element(window.ch0Rises.begin(), window.ch0Rises.end());
    for (const PulsePair& pair : window.ncal1Pairs) {
        const double relativeNs = pair.riseNs - referenceTimeNs;
        if (relativeNs < config.windowMinNs || relativeNs >= config.windowMaxNs) {
            continue;
        }
        if (pair.totNs < config.totMinNs || pair.totNs >= config.totMaxNs) {
            continue;
        }
        pulses.push_back(StoredPulse{static_cast<float>(relativeNs), static_cast<float>(pair.totNs)});
        if (passes_score_tot_gate(config, pulses.back().totNs)) {
            if ((scoreOrdinal % config.scoreStride) == 0) {
                const double globalTimeSeconds = window.hasGlobalTriggerSeconds
                    ? window.globalTriggerSeconds + pair.riseNs * 1.0e-9
                    : pair.riseNs * 1.0e-9;
                scorePulses.push_back(ScoredPulse{pulses.back().timeNs, globalTimeSeconds});
            }
            scoreOrdinal += 1;
        }
    }
}

CandidateMetrics evaluate_candidate(const Config& config,
                                    const std::vector<ScoredPulse>& scorePulses,
                                    double periodNs,
                                    int roundIndex) {
    CandidateMetrics metrics;
    metrics.periodNs = periodNs;
    metrics.roundIndex = roundIndex;
    if (scorePulses.empty()) {
        return metrics;
    }

    const int phaseBins = std::max(16, static_cast<int>(std::llround(periodNs / config.phaseBinWidthNs)));
    const double phaseBinWidthNs = periodNs / static_cast<double>(phaseBins);
    std::vector<std::uint64_t> phaseCounts(static_cast<std::size_t>(phaseBins), 0);
    double sumCos = 0.0;
    double sumSin = 0.0;

    for (const ScoredPulse& pulse : scorePulses) {
        const double timeNs = static_cast<double>(pulse.timeNs);
        const double phaseNs = positive_mod(timeNs, periodNs);
        const int phaseIndex = std::min(phaseBins - 1, static_cast<int>(std::floor(phaseNs / phaseBinWidthNs)));
        phaseCounts[static_cast<std::size_t>(phaseIndex)] += 1;
        const double angle = kTwoPi * phaseNs / periodNs;
        sumCos += std::cos(angle);
        sumSin += std::sin(angle);
    }

    const auto peakIt = std::max_element(phaseCounts.begin(), phaseCounts.end());
    const int peakIndex = static_cast<int>(std::distance(phaseCounts.begin(), peakIt));
    const double peakCenterNs = (static_cast<double>(peakIndex) + 0.5) * phaseBinWidthNs;
    metrics.peakCenterNs = peakCenterNs;
    metrics.peakHeight = peakIt == phaseCounts.end() ? 0 : *peakIt;
    metrics.coherence = std::sqrt(sumCos * sumCos + sumSin * sumSin) / static_cast<double>(scorePulses.size());

    double phaseOriginNs = peakCenterNs;
    for (int iteration = 0; iteration < 2; ++iteration) {
        double residualSum = 0.0;
        std::uint64_t selected = 0;
        for (const ScoredPulse& pulse : scorePulses) {
            const double phaseNs = positive_mod(static_cast<double>(pulse.timeNs), periodNs);
            const double residualNs = wrap_centered(phaseNs - phaseOriginNs, periodNs);
            if (std::abs(residualNs) <= config.peakWindowNs) {
                residualSum += residualNs;
                selected += 1;
            }
        }
        if (selected == 0) {
            return metrics;
        }
        phaseOriginNs = positive_mod(phaseOriginNs + residualSum / static_cast<double>(selected), periodNs);
    }

    double residualSum = 0.0;
    double residualSqSum = 0.0;
    double cycleSum = 0.0;
    double cycleSqSum = 0.0;
    double cycleResidualSum = 0.0;
    std::uint64_t selected = 0;
    for (const ScoredPulse& pulse : scorePulses) {
        const double timeNs = static_cast<double>(pulse.timeNs);
        const double phaseNs = positive_mod(timeNs, periodNs);
        const double residualNs = wrap_centered(phaseNs - phaseOriginNs, periodNs);
        if (std::abs(residualNs) > config.peakWindowNs) {
            continue;
        }
        const double cycleIndex = std::round((timeNs - phaseOriginNs) / periodNs);
        residualSum += residualNs;
        residualSqSum += residualNs * residualNs;
        cycleSum += cycleIndex;
        cycleSqSum += cycleIndex * cycleIndex;
        cycleResidualSum += cycleIndex * residualNs;
        selected += 1;
    }

    if (selected < config.minSelectedPulses) {
        return metrics;
    }

    const double selectedDouble = static_cast<double>(selected);
    const double meanResidualNs = residualSum / selectedDouble;
    const double varianceNs = std::max(0.0, residualSqSum / selectedDouble - meanResidualNs * meanResidualNs);
    metrics.phaseOriginNs = phaseOriginNs;
    metrics.selectedPulses = selected;
    metrics.selectedFraction = selectedDouble / static_cast<double>(scorePulses.size());
    metrics.meanResidualNs = meanResidualNs;
    metrics.sigmaNs = std::sqrt(varianceNs);

    const double denominator = selectedDouble * cycleSqSum - cycleSum * cycleSum;
    if (std::abs(denominator) > 1.0e-12) {
        metrics.driftSlopeNsPerCycle = (selectedDouble * cycleResidualSum - cycleSum * residualSum) / denominator;
        metrics.driftInterceptNs = (residualSum - metrics.driftSlopeNsPerCycle * cycleSum) / selectedDouble;
    }

    metrics.merit = metrics.selectedFraction / std::max(metrics.sigmaNs, 1.0e-12);
    metrics.valid = true;
    return metrics;
}

std::vector<CandidateMetrics> run_multistage_scan(const Config& config,
                                                  const std::vector<ScoredPulse>& scorePulses) {
    std::vector<CandidateMetrics> candidates;
    double scanMin = config.initialPeriodMinNs;
    double scanMax = config.initialPeriodMaxNs;
    double stepNs = config.initialStepNs;

    for (int roundIndex = 0; roundIndex < config.refineRounds; ++roundIndex) {
        CandidateMetrics roundBest;
        for (double periodNs = scanMin; periodNs <= scanMax + 0.5 * stepNs; periodNs += stepNs) {
            const CandidateMetrics metrics = evaluate_candidate(config, scorePulses, periodNs, roundIndex);
            candidates.push_back(metrics);
            if (metrics.valid && (!roundBest.valid || metrics.merit > roundBest.merit)) {
                roundBest = metrics;
            }
        }

        if (!roundBest.valid) {
            break;
        }

        scanMin = roundBest.periodNs - static_cast<double>(config.refineHalfSpanSteps) * stepNs;
        scanMax = roundBest.periodNs + static_cast<double>(config.refineHalfSpanSteps) * stepNs;
        stepNs /= config.refineFactor;
    }

    return candidates;
}

std::vector<CandidateMetrics> valid_candidates(const std::vector<CandidateMetrics>& candidates) {
    std::vector<CandidateMetrics> filtered;
    filtered.reserve(candidates.size());
    for (const auto& candidate : candidates) {
        if (candidate.valid) {
            filtered.push_back(candidate);
        }
    }
    return filtered;
}

Histogram1D build_best_phase_profile(const Config& config,
                                     const std::vector<ScoredPulse>& scorePulses,
                                     const CandidateMetrics& best) {
    const int phaseBins = std::max(16, static_cast<int>(std::llround(best.periodNs / config.phaseBinWidthNs)));
    Histogram1D histogram = make_histogram_1d(
        "custom_NcalRfPhaseProfile",
        "NCAL RF-phase profile from ch0-referenced times",
        phaseBins,
        0.0,
        best.periodNs,
        "RF phase (ns)"
    );
    for (const ScoredPulse& pulse : scorePulses) {
        const double phaseNs = positive_mod(static_cast<double>(pulse.timeNs) - best.phaseOriginNs, best.periodNs);
        fill_histogram(histogram, phaseNs);
    }
    return histogram;
}

Histogram2D build_best_phase_tot_histogram(const Config& config,
                                           const std::vector<StoredPulse>& pulses,
                                           const CandidateMetrics& best) {
    const int phaseBins = std::max(16, static_cast<int>(std::llround(best.periodNs / config.phaseBinWidthNs)));
    Histogram2D histogram = make_histogram_2d(
        "custom_NcalRfPhaseTot",
        "NCAL ToT vs recovered RF phase",
        phaseBins,
        0.0,
        best.periodNs,
        "RF phase (ns)",
        config.totBins,
        config.totMinNs,
        config.totMaxNs,
        "ToT (ns)"
    );
    for (const StoredPulse& pulse : pulses) {
        const double phaseNs = positive_mod(static_cast<double>(pulse.timeNs) - best.phaseOriginNs, best.periodNs);
        fill_histogram(histogram, phaseNs, static_cast<double>(pulse.totNs));
    }
    return histogram;
}

std::vector<CycleResidualStats> build_cycle_residuals(const Config& config,
                                                      const std::vector<ScoredPulse>& scorePulses,
                                                      const CandidateMetrics& best) {
    struct Accumulator {
        std::uint64_t count = 0;
        double sum = 0.0;
        double sumSq = 0.0;
    };

    std::map<std::int64_t, Accumulator> accumulators;
    for (const ScoredPulse& pulse : scorePulses) {
        const double timeNs = static_cast<double>(pulse.timeNs);
        const double phaseNs = positive_mod(timeNs - best.phaseOriginNs, best.periodNs);
        const double residualNs = wrap_centered(phaseNs, best.periodNs);
        if (std::abs(residualNs) > config.peakWindowNs) {
            continue;
        }
        const std::int64_t cycleIndex = static_cast<std::int64_t>(std::llround((timeNs - best.phaseOriginNs) / best.periodNs));
        auto& accumulator = accumulators[cycleIndex];
        accumulator.count += 1;
        accumulator.sum += residualNs;
        accumulator.sumSq += residualNs * residualNs;
    }

    std::vector<CycleResidualStats> rows;
    rows.reserve(accumulators.size());
    for (const auto& [cycleIndex, accumulator] : accumulators) {
        const double mean = accumulator.sum / static_cast<double>(accumulator.count);
        const double variance = std::max(0.0, accumulator.sumSq / static_cast<double>(accumulator.count) - mean * mean);
        rows.push_back(CycleResidualStats{cycleIndex, accumulator.count, mean, std::sqrt(variance)});
    }
    return rows;
}

std::vector<CycleResidualPoint> build_cycle_residual_points(const Config& config,
                                                            const std::vector<ScoredPulse>& scorePulses,
                                                            const CandidateMetrics& best) {
    std::vector<CycleResidualPoint> points;
    points.reserve(scorePulses.size());
    for (const ScoredPulse& pulse : scorePulses) {
        const double timeNs = static_cast<double>(pulse.timeNs);
        const double phaseNs = positive_mod(timeNs - best.phaseOriginNs, best.periodNs);
        const double residualNs = wrap_centered(phaseNs, best.periodNs);
        if (std::abs(residualNs) > config.peakWindowNs) {
            continue;
        }
        const std::int64_t cycleIndex = static_cast<std::int64_t>(std::llround((timeNs - best.phaseOriginNs) / best.periodNs));
        points.push_back(CycleResidualPoint{cycleIndex, residualNs, pulse.globalTimeSeconds});
    }
    return points;
}

void write_scan_file(const std::string& path,
                     const Config& config,
                     const std::vector<CandidateMetrics>& candidates) {
    std::ofstream output(path);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write scan file: " + path);
    }
    output << "# initial_period_range_ns=[" << format_double(config.initialPeriodMinNs) << ", "
           << format_double(config.initialPeriodMaxNs) << "]\n";
    output << "# initial_step_ns=" << format_double(config.initialStepNs) << '\n';
    output << "# refine_rounds=" << config.refineRounds << '\n';
    output << "# phase_bin_width_ns=" << format_double(config.phaseBinWidthNs) << '\n';
    output << "# peak_window_ns=" << format_double(config.peakWindowNs) << '\n';
    output << "# score_tot_gate_ns=[" << format_double(config.scoreTotMinNs) << ", "
           << format_double(config.scoreTotMaxNs) << ")\n";
    output << "# columns: round_index period_ns phase_origin_ns peak_center_ns peak_height selected_pulses selected_fraction sigma_ns mean_residual_ns coherence drift_slope_ns_per_cycle drift_intercept_ns merit valid\n";
    for (const CandidateMetrics& candidate : candidates) {
        output << candidate.roundIndex << ' '
               << format_double(candidate.periodNs) << ' '
               << format_double(candidate.phaseOriginNs) << ' '
               << format_double(candidate.peakCenterNs) << ' '
               << candidate.peakHeight << ' '
               << candidate.selectedPulses << ' '
               << format_double(candidate.selectedFraction) << ' '
               << format_double(candidate.sigmaNs) << ' '
               << format_double(candidate.meanResidualNs) << ' '
               << format_double(candidate.coherence) << ' '
               << format_double(candidate.driftSlopeNsPerCycle) << ' '
               << format_double(candidate.driftInterceptNs) << ' '
               << format_double(candidate.merit) << ' '
               << (candidate.valid ? "true" : "false") << '\n';
    }
}

void write_cycle_residual_file(const std::string& path,
                               const CandidateMetrics& best,
                               const std::vector<CycleResidualStats>& rows) {
    std::ofstream output(path);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write cycle residual file: " + path);
    }
    output << "# best_period_ns=" << format_double(best.periodNs) << '\n';
    output << "# phase_origin_ns=" << format_double(best.phaseOriginNs) << '\n';
    output << "# columns: cycle_index count mean_residual_ns rms_residual_ns\n";
    for (const CycleResidualStats& row : rows) {
        output << row.cycleIndex << ' '
               << row.count << ' '
               << format_double(row.meanResidualNs) << ' '
               << format_double(row.rmsResidualNs) << '\n';
    }
}

void write_cycle_residual_points_file(const std::string& path,
                                      const CandidateMetrics& best,
                                      const std::vector<CycleResidualPoint>& points) {
    std::ofstream output(path);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write cycle residual points file: " + path);
    }
    output << "# best_period_ns=" << format_double(best.periodNs) << '\n';
    output << "# phase_origin_ns=" << format_double(best.phaseOriginNs) << '\n';
    output << "# columns: cycle_index residual_ns global_time_seconds\n";
    for (const CycleResidualPoint& point : points) {
        output << point.cycleIndex << ' '
               << format_double(point.residualNs) << ' '
               << format_double(point.globalTimeSeconds) << '\n';
    }
}

void write_summary(const std::string& path,
                   const Config& config,
                   const AnalysisStats& stats,
                   std::uint64_t storedPulses,
                   std::uint64_t scorePulses,
                   std::uint64_t totalScoreEligiblePulses,
                   const CandidateMetrics& bestByMerit,
                   const CandidateMetrics& bestBySigma,
                   double elapsedSeconds) {
    std::ofstream output(path);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write summary file: " + path);
    }

    output << "input_file=" << config.inputPath << '\n';
    output << "window_min_ns=" << format_double(config.windowMinNs) << '\n';
    output << "window_max_ns=" << format_double(config.windowMaxNs) << '\n';
    output << "tot_gate_ns=[" << format_double(config.totMinNs) << ", " << format_double(config.totMaxNs) << ")\n";
    output << "score_tot_gate_ns=[" << format_double(config.scoreTotMinNs) << ", " << format_double(config.scoreTotMaxNs) << ")\n";
    output << "initial_period_range_ns=[" << format_double(config.initialPeriodMinNs) << ", " << format_double(config.initialPeriodMaxNs) << "]\n";
    output << "initial_step_ns=" << format_double(config.initialStepNs) << '\n';
    output << "refine_rounds=" << config.refineRounds << '\n';
    output << "refine_half_span_steps=" << config.refineHalfSpanSteps << '\n';
    output << "refine_factor=" << format_double(config.refineFactor) << '\n';
    output << "phase_bin_width_ns=" << format_double(config.phaseBinWidthNs) << '\n';
    output << "peak_window_ns=" << format_double(config.peakWindowNs) << '\n';
    output << "min_selected_pulses=" << config.minSelectedPulses << '\n';
    output << "score_stride=" << config.scoreStride << '\n';
    output << "parsed_headers=" << stats.parsedHeaders << '\n';
    output << "parsed_rows=" << stats.parsedRows << '\n';
    output << "finalized_windows=" << stats.finalizedWindows << '\n';
    output << "duplicate_restart_detected=" << (stats.duplicateRestartDetected ? "true" : "false") << '\n';
    output << "duplicate_restart_previous_trigger_seconds=" << format_double(stats.duplicateRestartPreviousTriggerSeconds) << '\n';
    output << "duplicate_restart_new_trigger_seconds=" << format_double(stats.duplicateRestartNewTriggerSeconds) << '\n';
    output << "duplicate_restart_headers_skipped=" << stats.duplicateRestartHeadersSkipped << '\n';
    output << "valid_ch0_rises=" << stats.validCh0Rises << '\n';
    output << "valid_ncal1_rises=" << stats.validNcal1Rises << '\n';
    output << "valid_ncal1_tot_pairs=" << stats.validNcal1Pairs << '\n';
    output << "windows_with_valid_ch0=" << stats.windowsWithValidCh0 << '\n';
    output << "windows_without_ch0=" << stats.windowsWithoutCh0 << '\n';
    output << "windows_with_multiple_ch0=" << stats.windowsWithMultipleCh0 << '\n';
    output << "stored_ncal_pulses=" << storedPulses << '\n';
    output << "score_eligible_ncal_pulses=" << totalScoreEligiblePulses << '\n';
    output << "score_ncal_pulses=" << scorePulses << '\n';
    output << "elapsed_seconds=" << format_double(elapsedSeconds) << "\n\n";

    const auto write_candidate = [&](const std::string& label, const CandidateMetrics& candidate) {
        output << label << '\n';
        output << "  valid=" << (candidate.valid ? "true" : "false") << '\n';
        output << "  period_ns=" << format_double(candidate.periodNs) << '\n';
        output << "  phase_origin_ns=" << format_double(candidate.phaseOriginNs) << '\n';
        output << "  peak_center_ns=" << format_double(candidate.peakCenterNs) << '\n';
        output << "  peak_height=" << candidate.peakHeight << '\n';
        output << "  selected_pulses=" << candidate.selectedPulses << '\n';
        output << "  selected_fraction=" << format_double(candidate.selectedFraction) << '\n';
        output << "  sigma_ns=" << format_double(candidate.sigmaNs) << '\n';
        output << "  mean_residual_ns=" << format_double(candidate.meanResidualNs) << '\n';
        output << "  coherence=" << format_double(candidate.coherence) << '\n';
        output << "  drift_slope_ns_per_cycle=" << format_double(candidate.driftSlopeNsPerCycle) << '\n';
        output << "  drift_intercept_ns=" << format_double(candidate.driftInterceptNs) << '\n';
        output << "  merit=" << format_double(candidate.merit) << "\n\n";
    };

    write_candidate("best_by_merit", bestByMerit);
    write_candidate("best_by_sigma", bestBySigma);

    output << "recipe\n";
    output << "  1. Build ch0-referenced NCAL rise times from the raw DOGMA stream.\n";
    output << "  2. Scan trial RF periods on the raw ch0-referenced times without using pre-folded histograms.\n";
    output << "  3. For each trial period, estimate the dominant phase origin from the phase peak.\n";
    output << "  4. Measure blur by the residual sigma inside a fixed peak window and track the selected fraction.\n";
    output << "  5. Check cycle-to-cycle drift by fitting residual vs cycle index.\n";
    output << "  6. Choose the best period from the scan and build the final NCAL phase-vs-ToT map.\n";
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
        std::vector<StoredPulse> pulses;
        std::vector<ScoredPulse> scorePulses;
        std::uint64_t scoreOrdinal = 0;
        WindowState currentWindow;

        int currentTdcOrdinal = -1;
        bool haveLastTrigger0 = false;
        double lastTrigger0Seconds = 0.0;
        bool haveSignalWindow = false;

        auto finalize_current_signal_window = [&]() {
            if (currentTdcOrdinal == kSignalOrdinal && haveSignalWindow) {
                finalize_window(currentWindow, config, stats, pulses, scorePulses, scoreOrdinal);
                currentWindow = WindowState{};
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
                    if (haveLastTrigger0) {
                        currentWindow.hasGlobalTriggerSeconds = true;
                        currentWindow.globalTriggerSeconds = lastTrigger0Seconds;
                    }
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
            if (channel != kCh0Channel && channel != kNcal1Channel) {
                continue;
            }
            if (!is_valid_hit_time_ns(timeNs)) {
                continue;
            }

            if (isRising == 1) {
                if (channel == kCh0Channel) {
                    currentWindow.ch0Rises.push_back(timeNs);
                    stats.validCh0Rises += 1;
                } else if (channel == kNcal1Channel) {
                    currentWindow.pendingNcal1Rises.push_back(timeNs);
                    stats.validNcal1Rises += 1;
                }
                continue;
            }

            if (channel == kNcal1Channel && !currentWindow.pendingNcal1Rises.empty()) {
                const double riseNs = currentWindow.pendingNcal1Rises.front();
                currentWindow.pendingNcal1Rises.pop_front();
                const double totNs = timeNs - riseNs;
                if ((totNs > 0.0) && std::isfinite(totNs)) {
                    currentWindow.ncal1Pairs.push_back(PulsePair{riseNs, totNs});
                    stats.validNcal1Pairs += 1;
                }
            }
        }

        finalize_current_signal_window();

        if (scorePulses.empty()) {
            throw std::runtime_error("No NCAL pulses survived the requested scoring gate");
        }

        const std::vector<CandidateMetrics> allCandidates = run_multistage_scan(config, scorePulses);
        const std::vector<CandidateMetrics> validCandidates = valid_candidates(allCandidates);
        if (validCandidates.empty()) {
            throw std::runtime_error("No valid RF-period candidates were found");
        }

        const auto bestByMeritIt = std::max_element(validCandidates.begin(), validCandidates.end(), [](const CandidateMetrics& left, const CandidateMetrics& right) {
            return left.merit < right.merit;
        });
        const auto bestBySigmaIt = std::min_element(validCandidates.begin(), validCandidates.end(), [](const CandidateMetrics& left, const CandidateMetrics& right) {
            return left.sigmaNs < right.sigmaNs;
        });
        const CandidateMetrics bestByMerit = *bestByMeritIt;
        const CandidateMetrics bestBySigma = *bestBySigmaIt;

        const Histogram1D bestPhaseProfile = build_best_phase_profile(config, scorePulses, bestByMerit);
        const Histogram2D bestPhaseTot = build_best_phase_tot_histogram(config, pulses, bestByMerit);
        const std::vector<CycleResidualStats> cycleResiduals = build_cycle_residuals(config, scorePulses, bestByMerit);
        const std::vector<CycleResidualPoint> cycleResidualPoints = build_cycle_residual_points(config, scorePulses, bestByMerit);

        const auto endedAt = std::chrono::steady_clock::now();
        const double elapsedSeconds = std::chrono::duration_cast<std::chrono::duration<double>>(endedAt - startedAt).count();

        write_scan_file(config.outputPrefix + "_scan.txt", config, allCandidates);
        write_histogram_1d(config.outputPrefix + "_best_phase_profile_hist.txt", bestPhaseProfile);
        write_histogram_2d(config.outputPrefix + "_best_phase_tot_hist.txt", bestPhaseTot);
        write_cycle_residual_file(config.outputPrefix + "_best_cycle_residuals.txt", bestByMerit, cycleResiduals);
        write_cycle_residual_points_file(config.outputPrefix + "_best_cycle_residual_points.txt", bestByMerit, cycleResidualPoints);
        write_summary(config.outputPrefix + "_summary.txt", config, stats, pulses.size(), scorePulses.size(), scoreOrdinal, bestByMerit, bestBySigma, elapsedSeconds);

        std::cout << "Stored NCAL pulses: " << pulses.size() << '\n';
        std::cout << "Score-eligible NCAL pulses: " << scoreOrdinal << '\n';
        std::cout << "Score NCAL pulses: " << scorePulses.size() << '\n';
        std::cout << "Best-by-merit period (ns): " << format_double(bestByMerit.periodNs) << '\n';
        std::cout << "Best-by-merit sigma (ns): " << format_double(bestByMerit.sigmaNs) << '\n';
        std::cout << "Best-by-merit drift slope (ns/cycle): " << format_double(bestByMerit.driftSlopeNsPerCycle) << '\n';
        std::cout << "Best-by-sigma period (ns): " << format_double(bestBySigma.periodNs) << '\n';
        std::cout << "Elapsed seconds: " << format_double(elapsedSeconds) << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}