#include <algorithm>
#include <array>
#include <charconv>
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
#include <random>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace {

constexpr double kSentinelValue = 2748779069440.0;
constexpr std::size_t kReservoirSize = 200000;
constexpr std::size_t kExampleLimit = 8;

struct ReservoirSampler {
    std::vector<double> values;
    std::uint64_t seen = 0;
    std::mt19937_64 rng{0xD06AULL};

    void add(double value) {
        ++seen;
        if (values.size() < kReservoirSize) {
            values.push_back(value);
            return;
        }

        std::uniform_int_distribution<std::uint64_t> dist(0, seen - 1);
        const auto slot = dist(rng);
        if (slot < values.size()) {
            values[static_cast<std::size_t>(slot)] = value;
        }
    }

    std::optional<double> percentile(double fraction) const {
        if (values.empty()) {
            return std::nullopt;
        }
        auto sorted = values;
        std::sort(sorted.begin(), sorted.end());
        const auto index = static_cast<std::size_t>(fraction * static_cast<double>(sorted.size() - 1));
        return sorted[index];
    }
};

struct ChannelStats {
    std::uint64_t rows = 0;
    std::uint64_t rise = 0;
    std::uint64_t fall = 0;
    std::uint64_t paired = 0;
    std::uint64_t duplicateRise = 0;
    std::uint64_t orphanRiseEndBlock = 0;
    std::uint64_t orphanFallWithoutRise = 0;
    std::uint64_t sentinelRows = 0;
    std::uint64_t negativeTot = 0;
    std::uint64_t zeroTot = 0;
    double totMin = std::numeric_limits<double>::infinity();
    double totMax = -std::numeric_limits<double>::infinity();
    double totSum = 0.0;
    ReservoirSampler sampler;
    std::map<int, std::uint64_t> firstColCounts;
    std::map<std::string, std::uint64_t> tdcCounts;
    std::vector<std::string> examples;
};

struct PendingRise {
    bool hasValue = false;
    double value = 0.0;
};

struct Config {
    std::string inputPath;
    std::vector<int> channels{1, 3};
    std::string anomalyOutputPath;
    std::uint64_t maxBlocks = 0;
};

bool starts_with(std::string_view text, std::string_view prefix) {
    return text.size() >= prefix.size() && text.substr(0, prefix.size()) == prefix;
}

std::vector<int> parse_channels(const std::string& text) {
    std::vector<int> channels;
    std::stringstream stream(text);
    std::string token;
    while (std::getline(stream, token, ',')) {
        if (token.empty()) {
            continue;
        }
        channels.push_back(std::stoi(token));
    }
    if (channels.empty()) {
        throw std::runtime_error("No channels were parsed from --channels");
    }
    return channels;
}

Config parse_args(int argc, char** argv) {
    Config config;
    for (int index = 1; index < argc; ++index) {
        const std::string arg = argv[index];
        if (arg == "--input" && index + 1 < argc) {
            config.inputPath = argv[++index];
        } else if (arg == "--channels" && index + 1 < argc) {
            config.channels = parse_channels(argv[++index]);
        } else if (arg == "--anomaly-output" && index + 1 < argc) {
            config.anomalyOutputPath = argv[++index];
        } else if (arg == "--max-blocks" && index + 1 < argc) {
            config.maxBlocks = std::stoull(argv[++index]);
        } else if (arg == "--help") {
            std::cout
                << "Usage: dogma_channel_scan --input <file> [--channels 1,3] [--anomaly-output path] [--max-blocks N]\n"
                << "Scans DOGMA ASCII .dat files for rise/fall pairing anomalies and TOT statistics.\n";
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown or incomplete argument: " + arg);
        }
    }

    if (config.inputPath.empty()) {
        throw std::runtime_error("Missing required argument: --input <file>");
    }
    return config;
}

std::string parse_tdc_id(const std::string& line) {
    const std::string key = "TDC TDC_";
    const auto start = line.find(key);
    if (start == std::string::npos) {
        return "UNKNOWN";
    }
    const auto id_start = start + 4;
    auto id_end = line.find(' ', id_start);
    if (id_end == std::string::npos) {
        id_end = line.size();
    }
    return line.substr(id_start, id_end - id_start);
}

int parse_size(const std::string& line) {
    const std::string key = " size ";
    const auto start = line.rfind(key);
    if (start == std::string::npos) {
        return -1;
    }
    const char* first = line.c_str() + start + key.size();
    char* end = nullptr;
    const long value = std::strtol(first, &end, 10);
    return end == first ? -1 : static_cast<int>(value);
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

    const long a = std::strtol(cursor, &end, 10);
    if (end == cursor || *end != ' ') {
        return false;
    }
    cursor = end + 1;

    const long b = std::strtol(cursor, &end, 10);
    if (end == cursor || *end != ' ') {
        return false;
    }
    cursor = end + 1;

    const long c = std::strtol(cursor, &end, 10);
    if (end == cursor || *end != ' ') {
        return false;
    }
    cursor = end + 1;

    const double d = std::strtod(cursor, &end);
    if (end == cursor || *end != '\0') {
        return false;
    }

    firstCol = static_cast<int>(a);
    channel = static_cast<int>(b);
    edge = static_cast<int>(c);
    value = d;
    return true;
}

void add_example(ChannelStats& stats, const std::string& text) {
    if (stats.examples.size() < kExampleLimit) {
        stats.examples.push_back(text);
    }
}

std::string format_double(double value) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(6) << value;
    return stream.str();
}

void update_tot(ChannelStats& stats, double tot) {
    stats.paired += 1;
    stats.totMin = std::min(stats.totMin, tot);
    stats.totMax = std::max(stats.totMax, tot);
    stats.totSum += tot;
    stats.sampler.add(tot);
    if (tot < 0.0) {
        stats.negativeTot += 1;
    } else if (tot == 0.0) {
        stats.zeroTot += 1;
    }
}

void log_anomaly(std::ofstream* output,
                 ChannelStats& stats,
                 const std::string& text) {
    add_example(stats, text);
    if (output != nullptr) {
        *output << text << '\n';
    }
}

void print_counter(const std::map<int, std::uint64_t>& counter) {
    bool first = true;
    std::cout << "{";
    for (const auto& [key, value] : counter) {
        if (!first) {
            std::cout << ", ";
        }
        first = false;
        std::cout << key << ": " << value;
    }
    std::cout << "}";
}

void print_counter(const std::map<std::string, std::uint64_t>& counter) {
    bool first = true;
    std::cout << "{";
    for (const auto& [key, value] : counter) {
        if (!first) {
            std::cout << ", ";
        }
        first = false;
        std::cout << key << ": " << value;
    }
    std::cout << "}";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Config config = parse_args(argc, argv);
        std::ifstream input(config.inputPath);
        if (!input.is_open()) {
            throw std::runtime_error("Unable to open input file: " + config.inputPath);
        }

        std::ofstream anomalyOutput;
        if (!config.anomalyOutputPath.empty()) {
            anomalyOutput.open(config.anomalyOutputPath);
            if (!anomalyOutput.is_open()) {
                throw std::runtime_error("Unable to open anomaly output file: " + config.anomalyOutputPath);
            }
            anomalyOutput << "# source_file=" << config.inputPath << "\n";
            anomalyOutput << "# channels=";
            for (std::size_t index = 0; index < config.channels.size(); ++index) {
                if (index > 0) {
                    anomalyOutput << ',';
                }
                anomalyOutput << config.channels[index];
            }
            anomalyOutput << "\n";
            anomalyOutput << "# format=event_index type tdc trigger channel details\n";
        }

        std::map<int, ChannelStats> stats;
        std::map<int, PendingRise> pending;
        for (const int channel : config.channels) {
            stats.emplace(channel, ChannelStats{});
            pending.emplace(channel, PendingRise{});
        }

        std::uint64_t anomalyIndex = 0;

        auto flush_block = [&](const std::string& tdcId, std::optional<double> trigger) {
            for (const int channel : config.channels) {
                auto& pendingRise = pending[channel];
                if (!pendingRise.hasValue) {
                    continue;
                }
                auto& channelStats = stats[channel];
                channelStats.orphanRiseEndBlock += 1;
                std::ostringstream example;
                example << ++anomalyIndex << " block_end_orphan_rise tdc=" << tdcId << " trigger=";
                if (trigger.has_value()) {
                    example << std::fixed << std::setprecision(6) << *trigger;
                } else {
                    example << "NA";
                }
                example << " ch=" << channel << " rise=" << std::fixed << std::setprecision(6)
                        << pendingRise.value;
                log_anomaly(config.anomalyOutputPath.empty() ? nullptr : &anomalyOutput, channelStats, example.str());
                pendingRise = PendingRise{};
            }
        };

        std::string currentTdcId = "UNKNOWN";
        std::optional<double> currentTrigger;
        std::uint64_t completedBlocks = 0;
        std::string line;
        const auto startedAt = std::chrono::steady_clock::now();

        while (std::getline(input, line)) {
            if (starts_with(line, "TDC ")) {
                if (completedBlocks > 0 || currentTrigger.has_value() || currentTdcId != "UNKNOWN") {
                    flush_block(currentTdcId, currentTrigger);
                    completedBlocks += 1;
                    if (config.maxBlocks > 0 && completedBlocks >= config.maxBlocks) {
                        break;
                    }
                }
                currentTdcId = parse_tdc_id(line);
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
            double value = 0.0;
            if (!parse_row(line, firstCol, channel, edge, value)) {
                continue;
            }

            auto statIt = stats.find(channel);
            if (statIt == stats.end()) {
                continue;
            }

            auto& channelStats = statIt->second;
            channelStats.rows += 1;
            channelStats.firstColCounts[firstCol] += 1;
            channelStats.tdcCounts[currentTdcId] += 1;
            if (value == kSentinelValue) {
                channelStats.sentinelRows += 1;
            }

            auto& pendingRise = pending[channel];
            if (edge == 1) {
                channelStats.rise += 1;
                if (pendingRise.hasValue) {
                    channelStats.duplicateRise += 1;
                    std::ostringstream example;
                    example << ++anomalyIndex << " duplicate_rise_before_fall tdc=" << currentTdcId << " trigger=";
                    if (currentTrigger.has_value()) {
                        example << std::fixed << std::setprecision(6) << *currentTrigger;
                    } else {
                        example << "NA";
                    }
                    example << " ch=" << channel << " prev_rise=" << std::fixed << std::setprecision(6)
                            << pendingRise.value << " new_rise=" << value;
                    log_anomaly(config.anomalyOutputPath.empty() ? nullptr : &anomalyOutput, channelStats, example.str());
                }
                pendingRise = PendingRise{true, value};
                continue;
            }

            channelStats.fall += 1;
            if (!pendingRise.hasValue) {
                channelStats.orphanFallWithoutRise += 1;
                std::ostringstream example;
                example << ++anomalyIndex << " fall_without_rise tdc=" << currentTdcId << " trigger=";
                if (currentTrigger.has_value()) {
                    example << std::fixed << std::setprecision(6) << *currentTrigger;
                } else {
                    example << "NA";
                }
                example << " ch=" << channel << " fall=" << std::fixed << std::setprecision(6) << value;
                log_anomaly(config.anomalyOutputPath.empty() ? nullptr : &anomalyOutput, channelStats, example.str());
                continue;
            }

            const double tot = value - pendingRise.value;
            if (tot < 0.0) {
                std::ostringstream example;
                example << ++anomalyIndex << " negative_tot tdc=" << currentTdcId << " trigger=";
                if (currentTrigger.has_value()) {
                    example << std::fixed << std::setprecision(6) << *currentTrigger;
                } else {
                    example << "NA";
                }
                example << " ch=" << channel << " rise=" << std::fixed << std::setprecision(6)
                        << pendingRise.value << " fall=" << value << " tot=" << tot;
                log_anomaly(config.anomalyOutputPath.empty() ? nullptr : &anomalyOutput, channelStats, example.str());
            }
            update_tot(channelStats, tot);
            pendingRise = PendingRise{};
        }

        if (config.maxBlocks == 0 || completedBlocks < config.maxBlocks) {
            flush_block(currentTdcId, currentTrigger);
            completedBlocks += 1;
        }

        const auto endedAt = std::chrono::steady_clock::now();
        const auto elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(endedAt - startedAt).count();

        std::cout << "File: " << config.inputPath << "\n";
        if (!config.anomalyOutputPath.empty()) {
            std::cout << "Anomaly log: " << config.anomalyOutputPath << "\n";
        }
        std::cout << "Completed TDC blocks: " << completedBlocks << "\n";
        std::cout << "Elapsed seconds: " << std::fixed << std::setprecision(3) << elapsed << "\n";

        for (const int channel : config.channels) {
            const auto& channelStats = stats[channel];
            std::cout << "\nChannel " << channel << "\n";
            std::cout << "  rows=" << channelStats.rows
                      << " rise=" << channelStats.rise
                      << " fall=" << channelStats.fall
                      << " paired=" << channelStats.paired << "\n";
            std::cout << "  sentinel_rows=" << channelStats.sentinelRows << "\n";
            std::cout << "  duplicate_rise=" << channelStats.duplicateRise << "\n";
            std::cout << "  orphan_rise_end_block=" << channelStats.orphanRiseEndBlock << "\n";
            std::cout << "  orphan_fall_without_rise=" << channelStats.orphanFallWithoutRise << "\n";
            std::cout << "  negative_tot=" << channelStats.negativeTot
                      << " zero_tot=" << channelStats.zeroTot << "\n";

            if (channelStats.paired > 0) {
                const auto p01 = channelStats.sampler.percentile(0.01).value();
                const auto p05 = channelStats.sampler.percentile(0.05).value();
                const auto p50 = channelStats.sampler.percentile(0.50).value();
                const auto p95 = channelStats.sampler.percentile(0.95).value();
                const auto p99 = channelStats.sampler.percentile(0.99).value();
                const double mean = channelStats.totSum / static_cast<double>(channelStats.paired);
                std::cout << "  TOT count=" << channelStats.paired
                          << " min=" << format_double(channelStats.totMin)
                          << " p01~=" << format_double(p01)
                          << " p05~=" << format_double(p05)
                          << " median~=" << format_double(p50)
                          << " p95~=" << format_double(p95)
                          << " p99~=" << format_double(p99)
                          << " max=" << format_double(channelStats.totMax)
                          << " mean=" << format_double(mean) << "\n";
            } else {
                std::cout << "  TOT count=0\n";
            }

            std::cout << "  first-column values=";
            print_counter(channelStats.firstColCounts);
            std::cout << "\n";
            std::cout << "  TDC ids=";
            print_counter(channelStats.tdcCounts);
            std::cout << "\n";
            for (const auto& example : channelStats.examples) {
                std::cout << "   - " << example << "\n";
            }
        }
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << "\n";
        return 1;
    }

    return 0;
}