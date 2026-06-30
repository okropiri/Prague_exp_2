#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
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
    int eventStride = 1;
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

struct EventState {
    std::array<std::array<std::array<double, 2>, 32>, 2> times{};
    bool hasAnyData = false;
};

struct AnalysisStats {
    std::uint64_t parsedHeaders = 0;
    std::uint64_t parsedRows = 0;
    std::uint64_t finalizedEvents = 0;
    std::uint64_t processedEvents = 0;
    std::uint64_t duplicateRestartHeadersSkipped = 0;
    bool duplicateRestartDetected = false;
    double duplicateRestartPreviousTriggerSeconds = 0.0;
    double duplicateRestartNewTriggerSeconds = 0.0;
};

constexpr int kNcal1Channel = 2;
constexpr int kSmallStilbeneChannel = 20;
constexpr int kLargeStilbeneChannel = 22;
constexpr int kRfChannel = 0;
constexpr int kTargetTdcIndex = 1;

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
        } else if (arg == "--event-stride" && index + 1 < argc) {
            config.eventStride = std::stoi(argv[++index]);
        } else if (arg == "--trigger-reset-threshold-seconds" && index + 1 < argc) {
            config.triggerResetThresholdSeconds = std::stod(argv[++index]);
        } else if (arg == "--help") {
            std::cout
                << "Usage: dogma_trgref_colleague_style --input <file> --output-prefix <prefix>\n"
                << "  [--event-stride 1] [--trigger-reset-threshold-seconds 1.0]\n";
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
    if (config.eventStride <= 0) {
        throw std::runtime_error("event-stride must be positive");
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
                              int yBins,
                              double yMin,
                              double yMax,
                              const std::string& xTitle,
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

    int xIndex = 0;
    if (xValue >= histogram.xMax) {
        xIndex = histogram.xBins + 1;
    } else if (xValue >= histogram.xMin) {
        const double fraction = (xValue - histogram.xMin) / (histogram.xMax - histogram.xMin);
        xIndex = 1 + std::min(histogram.xBins - 1, static_cast<int>(std::floor(fraction * histogram.xBins)));
    }

    int yIndex = 0;
    if (yValue >= histogram.yMax) {
        yIndex = histogram.yBins + 1;
    } else if (yValue >= histogram.yMin) {
        const double fraction = (yValue - histogram.yMin) / (histogram.yMax - histogram.yMin);
        yIndex = 1 + std::min(histogram.yBins - 1, static_cast<int>(std::floor(fraction * histogram.yBins)));
    }

    const auto flatIndex = static_cast<std::size_t>(yIndex) * static_cast<std::size_t>(histogram.xBins + 2)
                         + static_cast<std::size_t>(xIndex);
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
    output << "# columns: x_index y_index count\n";
    for (int yIndex = 0; yIndex < histogram.yBins + 2; ++yIndex) {
        const auto rowOffset = static_cast<std::size_t>(yIndex) * static_cast<std::size_t>(histogram.xBins + 2);
        for (int xIndex = 0; xIndex < histogram.xBins + 2; ++xIndex) {
            const auto count = histogram.storage[rowOffset + static_cast<std::size_t>(xIndex)];
            if (count == 0) {
                continue;
            }
            output << xIndex << ' ' << yIndex << ' ' << count << '\n';
        }
    }
}

void write_summary(const std::string& path,
                   const Config& config,
                   const AnalysisStats& stats,
                   const Histogram1D& ncal1,
                   const Histogram2D& profNcal1,
                   const Histogram2D& profLstilbene,
                   const Histogram2D& profSstilbene,
                   double elapsedSeconds) {
    std::ofstream output(path);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write summary file: " + path);
    }
    output << "input_file=" << config.inputPath << '\n';
    output << "event_stride=" << config.eventStride << '\n';
    output << "trigger_reset_threshold_seconds=" << format_double(config.triggerResetThresholdSeconds) << '\n';
    output << "parsed_headers=" << stats.parsedHeaders << '\n';
    output << "parsed_rows=" << stats.parsedRows << '\n';
    output << "finalized_events=" << stats.finalizedEvents << '\n';
    output << "processed_events=" << stats.processedEvents << '\n';
    output << "duplicate_restart_detected=" << (stats.duplicateRestartDetected ? "true" : "false") << '\n';
    output << "duplicate_restart_previous_trigger_seconds=" << format_double(stats.duplicateRestartPreviousTriggerSeconds) << '\n';
    output << "duplicate_restart_new_trigger_seconds=" << format_double(stats.duplicateRestartNewTriggerSeconds) << '\n';
    output << "duplicate_restart_headers_skipped=" << stats.duplicateRestartHeadersSkipped << '\n';
    output << "elapsed_seconds=" << format_double(elapsedSeconds) << "\n\n";

    output << "custom_Ncal1_entries=" << format_double(ncal1.entries) << '\n';
    output << "custom_ProfNcal1_entries=" << format_double(profNcal1.entries) << '\n';
    output << "custom_ProfLstilbene_entries=" << format_double(profLstilbene.entries) << '\n';
    output << "custom_ProfSstilbene_entries=" << format_double(profSstilbene.entries) << '\n';
}

void finalize_event(const Config& config,
                    AnalysisStats& stats,
                    EventState& event,
                    Histogram1D& hNcal1,
                    Histogram2D& hProfNcal1,
                    Histogram2D& hProfLstilbene,
                    Histogram2D& hProfSstilbene) {
    if (!event.hasAnyData) {
        return;
    }
    stats.finalizedEvents += 1;
    if (stats.finalizedEvents % static_cast<std::uint64_t>(config.eventStride) != 0) {
        event = EventState{};
        return;
    }

    stats.processedEvents += 1;

    const double ncal1Time = event.times[kTargetTdcIndex][kNcal1Channel][1];
    const double ncal1Fall = event.times[kTargetTdcIndex][kNcal1Channel][0];
    fill_histogram(hNcal1, ncal1Time);
    if (ncal1Time != 0.0 && ncal1Fall != 0.0) {
        fill_histogram(hProfNcal1, ncal1Time, ncal1Fall - ncal1Time);
    }

    const double sStilbeneTime = event.times[kTargetTdcIndex][kSmallStilbeneChannel][1];
    const double sStilbeneFall = event.times[kTargetTdcIndex][kSmallStilbeneChannel][0];
    if (sStilbeneTime != 0.0 && sStilbeneFall != 0.0) {
        fill_histogram(hProfSstilbene, sStilbeneTime, sStilbeneFall - sStilbeneTime);
    }

    const double lStilbeneTime = event.times[kTargetTdcIndex][kLargeStilbeneChannel][1];
    const double lStilbeneFall = event.times[kTargetTdcIndex][kLargeStilbeneChannel][0];
    if (lStilbeneTime != 0.0 && lStilbeneFall != 0.0) {
        fill_histogram(hProfLstilbene, lStilbeneTime, lStilbeneFall - lStilbeneTime);
    }

    event = EventState{};
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Config config = parse_args(argc, argv);
        std::ifstream input(config.inputPath);
        if (!input.is_open()) {
            throw std::runtime_error("Unable to open input file: " + config.inputPath);
        }

        Histogram1D hNcal1 = make_histogram_1d("custom_Ncal1", "Ncal1_TDC1", 1000, -800.0, 200.0, "ns");
        Histogram2D hProfNcal1 = make_histogram_2d("custom_ProfNcal1", "Ncal1_ToT_T", 12000, -6000.0, 6000.0, 300, 0.0, 300.0, "T (ns)", "TOT (ns)");
        Histogram2D hProfLstilbene = make_histogram_2d("custom_ProfLstilbene", "Lstilbene_ToT_T", 12000, -6000.0, 6000.0, 300, 0.0, 300.0, "T (ns)", "TOT (ns)");
        Histogram2D hProfSstilbene = make_histogram_2d("custom_ProfSstilbene", "Sstilbene_ToT_T", 12000, -6000.0, 6000.0, 300, 0.0, 300.0, "T (ns)", "TOT (ns)");
        AnalysisStats stats;
        EventState currentEvent;

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
                if (tdcOrdinal == 0) {
                    finalize_event(config, stats, currentEvent, hNcal1, hProfNcal1, hProfLstilbene, hProfSstilbene);
                }
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
            if (currentTdcOrdinal < 0 || currentTdcOrdinal >= 2) {
                continue;
            }
            if (channel < 0 || channel >= 32) {
                continue;
            }
            if (isRising != 0 && isRising != 1) {
                continue;
            }

            currentEvent.times[static_cast<std::size_t>(currentTdcOrdinal)][static_cast<std::size_t>(channel)][static_cast<std::size_t>(isRising)] = timeNs;
            currentEvent.hasAnyData = true;
        }

        finalize_event(config, stats, currentEvent, hNcal1, hProfNcal1, hProfLstilbene, hProfSstilbene);

        const auto endedAt = std::chrono::steady_clock::now();
        const double elapsedSeconds = std::chrono::duration_cast<std::chrono::duration<double>>(endedAt - startedAt).count();

        write_histogram_1d(config.outputPrefix + "_custom_Ncal1_hist.txt", hNcal1);
        write_histogram_2d(config.outputPrefix + "_custom_ProfNcal1_hist.txt", hProfNcal1);
        write_histogram_2d(config.outputPrefix + "_custom_ProfLstilbene_hist.txt", hProfLstilbene);
        write_histogram_2d(config.outputPrefix + "_custom_ProfSstilbene_hist.txt", hProfSstilbene);
        write_summary(config.outputPrefix + "_summary.txt", config, stats, hNcal1, hProfNcal1, hProfLstilbene, hProfSstilbene, elapsedSeconds);

        std::cout << "Input file: " << config.inputPath << '\n';
        std::cout << "Output prefix: " << config.outputPrefix << '\n';
        std::cout << "Finalized events: " << stats.finalizedEvents << '\n';
        std::cout << "Processed events: " << stats.processedEvents << '\n';
        std::cout << "Duplicate restart detected: " << (stats.duplicateRestartDetected ? "yes" : "no") << '\n';
        if (stats.duplicateRestartDetected) {
            std::cout << "Duplicate restart previous trigger [s]: " << format_double(stats.duplicateRestartPreviousTriggerSeconds) << '\n';
            std::cout << "Duplicate restart new trigger [s]: " << format_double(stats.duplicateRestartNewTriggerSeconds) << '\n';
        }
        std::cout << "custom_Ncal1 entries: " << format_double(hNcal1.entries) << '\n';
        std::cout << "custom_ProfNcal1 entries: " << format_double(hProfNcal1.entries) << '\n';
        std::cout << "custom_ProfLstilbene entries: " << format_double(hProfLstilbene.entries) << '\n';
        std::cout << "custom_ProfSstilbene entries: " << format_double(hProfSstilbene.entries) << '\n';
        std::cout << "Elapsed seconds: " << format_double(elapsedSeconds) << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}