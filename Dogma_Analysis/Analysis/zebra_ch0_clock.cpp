#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
constexpr int channels = 32;

struct Pulse {
    std::uint64_t window;
    int ordinal, channel, riseIndex, fallIndex;
    double trigger, rise, fall, tot;
};

struct Reference {
    Pulse pulse;
    int candidates = 0;
    std::int64_t tick = 0;
    bool triggerSeen = false, riseSeen = false, fallSeen = false;
};

using Options = std::map<std::string, std::string>;

double option(const Options& args, const std::string& name, double fallback) {
    const auto found = args.find(name);
    const double value = found == args.end() ? fallback : std::stod(found->second);
    if (!std::isfinite(value)) throw std::runtime_error("Nonfinite option " + name);
    return value;
}

bool allow_gaps(const Options& args) {
    const auto found = args.find("--allow-reference-gaps");
    if (found == args.end() || found->second == "0") return false;
    if (found->second == "1") return true;
    throw std::runtime_error("--allow-reference-gaps must be 0 or 1");
}

struct ReferenceGate {
    double riseMin, riseMax, totMin, totMax;
    explicit ReferenceGate(const Options& args)
        : riseMin(option(args, "--ch0-valid-rise-min-ns", -392)),
          riseMax(option(args, "--ch0-valid-rise-max-ns", -374)),
          totMin(option(args, "--ch0-valid-tot-min-ns", 16.5)),
          totMax(option(args, "--ch0-valid-tot-max-ns", 19.5)) {
        if (riseMin > riseMax || totMin <= 0 || totMin > totMax)
            throw std::runtime_error("Invalid ch0 reference gates");
    }
};

std::ifstream open_input(const fs::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot read " + path.string());
    return input;
}

std::ofstream open_output(const fs::path& path) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Cannot write " + path.string());
    output << std::setprecision(18);
    return output;
}

bool parse_pulse(const std::string& line, Pulse& p) {
    if (line.empty() || line[0] == '#' || line.rfind("window_index", 0) == 0) return false;
    std::istringstream stream(line);
    std::string extra;
    if (!(stream >> p.window >> p.ordinal >> p.trigger >> p.channel >> p.riseIndex
          >> p.fallIndex >> p.rise >> p.fall >> p.tot) || (stream >> extra)
        || p.window == 0 || p.channel < 0 || p.channel >= channels
        || !std::isfinite(p.trigger) || !std::isfinite(p.rise)
        || !std::isfinite(p.fall) || !std::isfinite(p.tot) || p.tot <= 0) {
        throw std::runtime_error("Invalid cleaned pulse row: " + line);
    }
    return true;
}

bool valid_reference(const Pulse& p, const ReferenceGate& gate) {
    return p.channel == 0 && p.rise >= gate.riseMin && p.rise <= gate.riseMax
        && p.tot >= gate.totMin && p.tot <= gate.totMax;
}

unsigned long hardware_number(std::string label) {
    if (label.rfind("TDC_", 0) == 0) label.erase(0, 4);
    std::size_t parsed = 0;
    const auto value = std::stoul(label, &parsed, 16);
    if (parsed != label.size()) throw std::runtime_error("Invalid hardware identity " + label);
    return value;
}

void extract(const Options& args) {
    const int ordinal = std::stoi(args.at("--ordinal"));
    const long double quantum = std::stold(args.at("--coarse-quantum-ns"));
    if (ordinal < 0 || !(quantum > 0 && quantum < 100)) throw std::runtime_error("Invalid TDC ordinal/coarse quantum");
    const ReferenceGate gate(args);
    const bool gaps = allow_gaps(args);
    auto input = open_input(args.at("--input"));
    std::deque<Reference> refs;
    struct MissingWindow {
        std::uint64_t window, pulses, ch0;
        std::string reason;
    };
    std::vector<MissingWindow> missing;
    std::array<std::uint64_t, channels> counts{};
    std::uint64_t current = 0, windowPulses = 0, windowCh0 = 0, represented = 0;
    Reference ref{};
    auto finish = [&]() {
        if (!current) return;
        ++represented;
        if (!ref.candidates) {
            if (!gaps) throw std::runtime_error("No guarded ch0 in window " + std::to_string(current));
            missing.push_back({current, windowPulses, windowCh0,
                               windowCh0 ? "no_passing_ch0" : "no_accepted_ch0"});
        } else {
            refs.push_back(ref);
        }
    };
    std::string line;
    Pulse p{};
    while (std::getline(input, line)) {
        if (!parse_pulse(line, p)) continue;
        if (p.ordinal != ordinal || p.riseIndex != ordinal + 1 || p.fallIndex != ordinal + 1)
            throw std::runtime_error("Unexpected TDC identity");
        if (p.window != current) {
            if (p.window <= current || (!gaps && p.window != current + 1))
                throw std::runtime_error("Non-contiguous or reordered cleaned windows");
            finish();
            for (auto empty = current + 1; empty < p.window; ++empty)
                missing.push_back({empty, 0, 0, "empty_cleaned_window"});
            current = p.window;
            ref = Reference{};
            windowPulses = windowCh0 = 0;
        }
        counts[p.channel]++;
        ++windowPulses;
        if (p.channel == 0) ++windowCh0;
        if (valid_reference(p, gate)) {
            if (!ref.candidates) ref.pulse = p;
            ref.candidates++;
        }
    }
    finish();
    if (refs.size() < 32) throw std::runtime_error("Insufficient ch0 references");

    auto decoded = open_input(args.at("--decoded"));
    std::uint64_t window = 0, board0Headers = 0;
    std::size_t referenceIndex = 0;
    bool target = false;
    Reference* selected = nullptr;
    long double maxLatticeError = 0;
    while (std::getline(decoded, line)) {
        if (line.rfind("TDC ", 0) == 0) {
            target = false;
            selected = nullptr;
            if (line.rfind("TDC 0 of ", 0) == 0) ++board0Headers;
            const std::string prefix = "TDC " + std::to_string(ordinal) + " of ";
            target = line.rfind(prefix, 0) == 0;
            if (target) {
                ++window;
                const auto hardware = args.find("--hardware-id");
                if (hardware != args.end()) {
                    const auto position = line.rfind(" TDC ");
                    if (position == std::string::npos) throw std::runtime_error("Missing decoded hardware identity");
                    std::istringstream label(line.substr(position + 5));
                    std::string value;
                    label >> value;
                    if (hardware_number(value) != hardware_number(hardware->second))
                        throw std::runtime_error("Decoded hardware identity mismatch");
                }
                if (referenceIndex < refs.size()) {
                    if (refs[referenceIndex].pulse.window < window)
                        throw std::runtime_error("Decoded/cleaned reference sequence mismatch");
                    if (refs[referenceIndex].pulse.window == window)
                        selected = &refs[referenceIndex++];
                }
            }
            continue;
        }
        if (!target || !selected) continue;
        Reference& r = *selected;
        if (line.rfind("GlobalTriggerTime ", 0) == 0) {
            if (r.triggerSeen) throw std::runtime_error("Duplicate trigger anchor");
            const std::string stamp = line.substr(18);
            const auto decimal = stamp.find('.');
            if (decimal == std::string::npos || stamp.size() - decimal - 1 < 9)
                throw std::runtime_error("Original >=9-decimal trigger metadata required");
            const long double ns = std::stold(stamp) * 1e9L;
            r.tick = std::llround(ns / quantum);
            const long double latticeError = std::abs(ns - r.tick * quantum);
            maxLatticeError = std::max(maxLatticeError, latticeError);
            if (latticeError > 0.65L)
                throw std::runtime_error("Trigger anchor is incompatible with decoder coarse lattice");
            if (std::abs(ns / 1e9L - r.pulse.trigger) > 0.501e-6L)
                throw std::runtime_error("Original trigger metadata does not match cleaned input");
            r.triggerSeen = true;
            continue;
        }
        if (line.rfind(std::to_string(ordinal + 1) + " 0 ", 0) != 0) continue;
        int rowOrdinal, channel, rising;
        double time;
        std::istringstream edge(line);
        if (!(edge >> rowOrdinal >> channel >> rising >> time)
            || rowOrdinal != ordinal + 1 || channel != 0)
            throw std::runtime_error("Reference edge identity mismatch");
        if (rising == 1 && std::abs(time - r.pulse.rise) <= 0.0015) {
            if (r.riseSeen) throw std::runtime_error("Duplicate reference rise");
            r.riseSeen = true;
        } else if (rising == 0 && std::abs(time - r.pulse.fall) <= 0.0015) {
            if (r.fallSeen) throw std::runtime_error("Duplicate reference fall");
            r.fallSeen = true;
        }
    }
    if (referenceIndex != refs.size() || window < current || (!gaps && window != current))
        throw std::runtime_error("Decoded event count mismatch");
    const auto expected = args.find("--expected-windows");
    if (expected != args.end() && window != std::stoull(expected->second))
        throw std::runtime_error("Decoded selected-board count differs from cleaning");
    for (auto empty = current + 1; empty <= window; ++empty)
        missing.push_back({empty, 0, 0, "empty_cleaned_window"});
    auto output = open_output(args.at("--output"));
    output << "# method=original_coarse_counter_plus_cleaned_fine_ch0\n"
           << "# input_pulse_table=" << args.at("--input") << '\n'
           << "# original_decoded_metadata=" << args.at("--decoded") << '\n'
           << "# coarse_quantum_ns=" << quantum << '\n'
           << "# max_trigger_lattice_rounding_error_ns=" << maxLatticeError << '\n'
           << "# decoded_selected_board_windows=" << window << '\n'
           << "# decoded_board0_headers=" << board0Headers << '\n'
           << "# represented_cleaned_windows=" << represented << '\n'
           << "# reference_gap_windows=" << missing.size() << '\n'
           << "# empty_cleaned_windows=" << window - represented << '\n'
           << "# input_channel_counts=";
    for (int ch = 0; ch < channels; ++ch) output << (ch ? "," : "") << counts[ch];
    output << "\n# columns=window_index coarse_tick rise_ns tot_ns elapsed_ns valid_candidates\n";
    for (const auto& r : refs) {
        if (!r.triggerSeen || !r.riseSeen || !r.fallSeen)
            throw std::runtime_error("Unverified reference in window " + std::to_string(r.pulse.window));
        // Subtract integer counters before converting to avoid large absolute-time cancellation.
        const long double elapsed = (r.tick - refs.front().tick) * quantum
            + static_cast<long double>(r.pulse.rise) - refs.front().pulse.rise;
        output << r.pulse.window << '\t' << r.tick << '\t' << r.pulse.rise << '\t'
               << r.pulse.tot << '\t' << elapsed << '\t' << r.candidates << '\n';
    }
    const fs::path gapPath = fs::path(args.at("--output")).parent_path() / "reference_gap_windows.tsv";
    auto gapOutput = open_output(gapPath);
    gapOutput << "window_index\treason\tcleaned_pulses\taccepted_ch0_pulses\n";
    for (const auto& gap : missing)
        gapOutput << gap.window << '\t' << gap.reason << '\t' << gap.pulses << '\t' << gap.ch0 << '\n';
    std::cout << "Verified " << refs.size() << " original/cleaned ch0 references\n";
}

using Storage = std::array<std::vector<std::uint64_t>, channels>;

struct Histograms {
    int xBins, yBins;
    double xMax, yMin, yMax;
    Storage data;

    Histograms(int nx, int ny, double xmax, double ymin, double ymax, bool detailed,
               const std::set<int>& selected)
        : xBins(nx), yBins(ny), xMax(xmax), yMin(ymin), yMax(ymax) {
        for (int ch = 0; ch < channels; ++ch)
            if (!detailed || selected.count(ch))
                data[ch].resize(static_cast<std::size_t>(nx) * ny);
    }

    bool fill(int ch, double x, double y = 0) {
        if (data[ch].empty()) return true;
        if (x < 0 || x >= xMax || y < yMin || y >= yMax) return false;
        const int xb = static_cast<int>(x / xMax * xBins);
        const int yb = static_cast<int>((y - yMin) / (yMax - yMin) * yBins);
        if (xb < 0 || xb >= xBins || yb < 0 || yb >= yBins)
            throw std::runtime_error("Histogram bin outside range");
        data[ch][static_cast<std::size_t>(yb) * xBins + xb]++;
        return true;
    }

    void write(const fs::path& path, bool profile, double period, const std::string& runKey,
               const std::string& periodModel) const {
        auto out = open_output(path);
        out << "# run_key=" << runKey << '\n'
            << "# x_min_ns=0\n# x_max_ns=" << xMax << "\n# x_bins=" << xBins
            << "\n# bins=" << xBins << "\n# channel_count=32\n# y_min_ns=" << yMin
            << "\n# y_max_ns=" << yMax << "\n# y_bins=" << yBins
            << "\n# deduced_period_ns=" << period
            << "\n# phase_origin_ns=0\n# period_model=" << periodModel
            << "\n# phase_axis=local_cycles_scaled_by_board_reference_period\n"
            << "# stored_channels=";
        bool first = true;
        for (int ch = 0; ch < channels; ++ch) if (!data[ch].empty()) {
            out << (first ? "" : ",") << ch;
            first = false;
        }
        out << '\n';
        if (profile) {
            out << "# columns=bin center ch0..ch31 total\n";
            for (int x = 0; x < xBins; ++x) {
                out << x + 1 << ' ' << (x + 0.5) * xMax / xBins;
                std::uint64_t total = 0;
                for (int ch = 0; ch < channels; ++ch) {
                    out << ' ' << data[ch][x];
                    total += data[ch][x];
                }
                out << ' ' << total << '\n';
            }
        } else {
            out << "# columns=channel x_bin y_bin count\n";
            for (int ch = 0; ch < channels; ++ch)
                for (std::size_t index = 0; index < data[ch].size(); ++index)
                    if (data[ch][index])
                        out << ch << ' ' << index % xBins << ' ' << index / xBins
                            << ' ' << data[ch][index] << '\n';
        }
    }
};

double positive_mod(double x, double period) {
    double value = std::fmod(x, period);
    if (value < 0) value += period;
    return value >= period ? 0 : value;
}

void fold(const Options& args) {
    const fs::path out = args.at("--output");
    if (fs::exists(out)) throw std::runtime_error("Refusing to replace existing folded output");
    fs::create_directory(out);
    const ReferenceGate gate(args);
    const bool gaps = allow_gaps(args);
    const double periodMin = option(args, "--period-min-ns", 35);
    const double allCh0Option = option(args, "--include-all-ch0-hits", 0);
    if (allCh0Option != 0 && allCh0Option != 1)
        throw std::runtime_error("--include-all-ch0-hits must be 0 or 1");
    const bool includeAllCh0 = allCh0Option == 1;
    const auto modelOption = args.find("--period-model");
    const std::string periodModel = modelOption == args.end()
        ? "independent_ch0_neighbor_interval_average" : modelOption->second;
    if (periodModel != "independent_ch0_neighbor_interval_average" && periodModel != "shared_burst_rf_period"
        && periodModel != "shared_ch2_detector_period")
        throw std::runtime_error("Unknown RF period model");
    const double periodMax = option(args, "--period-max-ns", 45);
    if (periodMin <= 0 || periodMin >= periodMax)
        throw std::runtime_error("Invalid clock period guard");
    std::set<std::uint64_t> unresolved;
    const auto unresolvedFile = args.find("--unresolved-clock-windows");
    if (unresolvedFile != args.end()) {
        if (!gaps) throw std::runtime_error("Unresolved clock windows require explicit gap handling");
        auto file = open_input(unresolvedFile->second);
        std::string row;
        while (std::getline(file, row)) {
            if (row.empty() || row[0] == '#') continue;
            std::istringstream values(row);
            std::uint64_t id;
            if (!(values >> id) || id == 0 || !unresolved.insert(id).second)
                throw std::runtime_error("Invalid/duplicate unresolved clock window");
        }
    }
    auto periodsFile = open_input(args.at("--periods"));
    std::vector<double> periods(1, 0);
    std::string line;
    std::uint64_t window;
    double period;
    std::uint64_t previous = 0, periodCount = 0;
    long double sum = 0;
    while (std::getline(periodsFile, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream row(line);
        if (!(row >> window >> period) || window <= previous || (!gaps && window != previous + 1)
            || !std::isfinite(period) || period < periodMin || period > periodMax || unresolved.count(window))
            throw std::runtime_error("Invalid per-window RF period");
        periods.resize(window + 1, 0);
        periods[window] = period;
        previous = window;
        ++periodCount;
        sum += period;
    }
    if (!periodCount) throw std::runtime_error("No calibrated periods");
    if (periodModel == "shared_burst_rf_period" || periodModel == "shared_ch2_detector_period") {
        double shared = 0;
        for (const auto value : periods) if (value > 0) {
            if (shared > 0 && value != shared) throw std::runtime_error("Shared RF-period model contains varying periods");
            shared = value;
        }
    }
    const double referencePeriod = sum / periodCount;
    const double phaseWidth = option(args, "--phase-bin-width-ns", 0.25);
    const double totMin = option(args, "--tot-min-ns", 0), totMax = option(args, "--tot-max-ns", 80);
    const double timeMin = option(args, "--ch0ref-min-ns", -11000), timeMax = option(args, "--ch0ref-max-ns", 2000);
    const double triggerMin = option(args, "--trigger-ref-min-ns", -11000), triggerMax = option(args, "--trigger-ref-max-ns", 2000);
    const int totBins = static_cast<int>(option(args, "--tot-bins", 800));
    const int timeBins = static_cast<int>(option(args, "--ch0ref-bins", 2600));
    const int triggerBins = static_cast<int>(option(args, "--trigger-ref-bins", 2600));
    if (phaseWidth <= 0 || totMin >= totMax || timeMin >= timeMax || triggerMin >= triggerMax
        || totBins < 1 || timeBins < 1 || triggerBins < 1)
        throw std::runtime_error("Invalid folding histogram bounds");
    std::set<int> detailed{0, 1, 2, 3, 6};
    const auto chosen = args.find("--histogram-channels");
    if (chosen != args.end()) {
        detailed.clear();
        std::istringstream values(chosen->second);
        std::string value;
        while (std::getline(values, value, ',')) {
            const int ch = std::stoi(value);
            if (ch < 0 || ch >= channels || !detailed.insert(ch).second)
                throw std::runtime_error("Invalid/duplicate detailed channel");
        }
        if (detailed.empty()) throw std::runtime_error("No detailed channels requested");
    }
    const int bins = static_cast<int>(std::ceil(referencePeriod / phaseWidth));
    Histograms profile(bins, 1, referencePeriod, 0, 1, false, detailed);
    Histograms profile3(3 * bins, 1, 3 * referencePeriod, 0, 1, false, detailed);
    Histograms tot(bins, totBins, referencePeriod, totMin, totMax, true, detailed);
    Histograms tot3(3 * bins, totBins, 3 * referencePeriod, totMin, totMax, true, detailed);
    Histograms time(bins, timeBins, referencePeriod, timeMin, timeMax, true, detailed);
    Histograms time3(3 * bins, timeBins, 3 * referencePeriod, timeMin, timeMax, true, detailed);
    Histograms trigger(bins, triggerBins, referencePeriod, triggerMin, triggerMax, true, detailed);
    Histograms trigger3(3 * bins, triggerBins, 3 * referencePeriod, triggerMin, triggerMax, true, detailed);
    std::array<std::uint64_t, channels> inputCounts{}, accepted{}, overflow{}, overflow80{}, withoutReference{}, nonreferenceCh0{}, withoutClock{};
    auto input = open_input(args.at("--input"));
    std::vector<Pulse> hits;
    std::uint64_t current = 0;
    std::uint64_t represented = 0, calibrated = 0, missing = 0, unresolvedSeen = 0;
    auto finish = [&]() {
        if (!current) return;
        ++represented;
        const Pulse* ref = nullptr;
        for (const auto& hit : hits) if (valid_reference(hit, gate) && !ref) ref = &hit;
        const double localPeriod = current < periods.size() ? periods[current] : 0;
        if (!ref) {
            if (unresolved.count(current)) throw std::runtime_error("Unresolved clock window has no reference");
            if (!gaps) throw std::runtime_error("Missing guarded ch0 while folding");
            if (localPeriod != 0) throw std::runtime_error("Period supplied for a window without a reference");
            ++missing;
            for (const auto& hit : hits) ++withoutReference[hit.channel];
            return;
        }
        if (localPeriod == 0) {
            if (!unresolved.count(current)) throw std::runtime_error("Reference-bearing window has no calibrated period");
            ++unresolvedSeen;
            for (const auto& hit : hits) ++withoutClock[hit.channel];
            return;
        }
        ++calibrated;
        for (const auto& hit : hits) {
            if (!includeAllCh0 && hit.channel == 0 && &hit != ref) {
                ++nonreferenceCh0[0];
                continue;
            }
            const double relative = hit.rise - ref->rise;
            const double cycles = relative / localPeriod;
            const double phase = positive_mod(cycles, 1) * referencePeriod;
            const double phase3 = positive_mod(cycles, 3) * referencePeriod;
            if (!profile.fill(hit.channel, phase) || !profile3.fill(hit.channel, phase3))
                throw std::runtime_error("Unexpected phase overflow");
            accepted[hit.channel]++;
            if (hit.tot < totMin || hit.tot >= totMax) overflow[hit.channel]++;
            if (hit.tot < 0 || hit.tot >= 80) overflow80[hit.channel]++;
            tot.fill(hit.channel, phase, hit.tot);
            tot3.fill(hit.channel, phase3, hit.tot);
            if (!time.fill(hit.channel, phase, relative) || !time3.fill(hit.channel, phase3, relative)
                || !trigger.fill(hit.channel, phase, hit.rise)
                || !trigger3.fill(hit.channel, phase3, hit.rise))
                throw std::runtime_error("Timing axes do not contain all detailed pulses");
        }
    };
    Pulse p{};
    while (std::getline(input, line)) {
        if (!parse_pulse(line, p)) continue;
        if (p.window != current) {
            if (p.window <= current || (!gaps && p.window != current + 1))
                throw std::runtime_error("Non-contiguous or reordered folding windows");
            finish();
            current = p.window;
            hits.clear();
        }
        const auto ordinal = args.find("--ordinal");
        if (ordinal != args.end()) {
            const int expectedOrdinal = std::stoi(ordinal->second);
            if (p.ordinal != expectedOrdinal || p.riseIndex != expectedOrdinal + 1 || p.fallIndex != expectedOrdinal + 1)
                throw std::runtime_error("Unexpected folding TDC identity");
        }
        inputCounts[p.channel]++;
        hits.push_back(p);
    }
    finish();
    if (calibrated != periodCount || unresolvedSeen != unresolved.size() || (!gaps && current != previous))
        throw std::runtime_error("Period/input event count mismatch");
    const auto expectedWindows = args.find("--expected-windows");
    const std::uint64_t totalWindows = expectedWindows == args.end() ? current : std::stoull(expectedWindows->second);
    if (current > totalWindows || previous > totalWindows)
        throw std::runtime_error("Folding window exceeds original selected-board count");
    const std::string key = args.at("--run-key");
    profile.write(out / "folded_matrix.tsv", true, referencePeriod, key, periodModel);
    profile3.write(out / "folded3_matrix.tsv", true, referencePeriod, key, periodModel);
    tot.write(out / "phase_tot.tsv", false, referencePeriod, key, periodModel);
    tot3.write(out / "phase3_tot.tsv", false, referencePeriod, key, periodModel);
    time.write(out / "phase_time.tsv", false, referencePeriod, key, periodModel);
    time3.write(out / "phase3_time.tsv", false, referencePeriod, key, periodModel);
    trigger.write(out / "phase_trigger.tsv", false, referencePeriod, key, periodModel);
    trigger3.write(out / "phase3_trigger.tsv", false, referencePeriod, key, periodModel);
    auto stats = open_output(out / "channel_counts.tsv");
    stats << "channel\tinput_pulses\tfolded_pulses\ttot_outside_0_80"
          << "\ttot_outside_display\twithout_reference_pulses\tnonreference_ch0_pulses\twithout_clock_pulses\n";
    for (int ch = 0; ch < channels; ++ch)
        stats << ch << '\t' << inputCounts[ch] << '\t' << accepted[ch] << '\t'
              << overflow80[ch] << '\t'
              << overflow[ch] << '\t' << withoutReference[ch] << '\t' << nonreferenceCh0[ch] << '\t' << withoutClock[ch] << '\n';
    auto summary = open_output(out / "fold_summary.txt");
    summary << "represented_windows=" << represented << "\nreference_windows=" << calibrated
            << "\nrepresented_windows_without_reference=" << missing
            << "\nunresolved_clock_reference_windows=" << unresolvedSeen
            << "\nempty_cleaned_windows=" << totalWindows - represented
            << "\nmean_period_ns=" << referencePeriod << '\n';
    std::cout << "Folded " << current << " windows, mean period " << referencePeriod << " ns\n";
}

int main(int argc, char** argv) {
    try {
        if (argc < 2 || (argc - 2) % 2) throw std::runtime_error("Usage: zebra_ch0_clock extract|fold --key value ...");
        Options args;
        for (int i = 2; i < argc; i += 2) {
            if (!args.emplace(argv[i], argv[i + 1]).second)
                throw std::runtime_error("Duplicate option " + std::string(argv[i]));
        }
        const std::string mode = argv[1];
        std::set<std::string> allowed{"--input", "--output", "--ch0-valid-rise-min-ns",
            "--ch0-valid-rise-max-ns", "--ch0-valid-tot-min-ns", "--ch0-valid-tot-max-ns",
            "--allow-reference-gaps", "--expected-windows", "--ordinal"};
        if (mode == "extract") {
            allowed.insert({"--decoded", "--coarse-quantum-ns", "--hardware-id"});
        } else if (mode == "fold") {
            allowed.insert({"--periods", "--run-key", "--phase-bin-width-ns", "--tot-min-ns",
                "--tot-max-ns", "--tot-bins", "--ch0ref-min-ns", "--ch0ref-max-ns",
                "--ch0ref-bins", "--trigger-ref-min-ns", "--trigger-ref-max-ns",
                "--trigger-ref-bins", "--histogram-channels", "--period-min-ns", "--period-max-ns",
                "--unresolved-clock-windows", "--include-all-ch0-hits", "--period-model"});
        } else {
            throw std::runtime_error("Unknown mode " + mode);
        }
        for (const auto& entry : args)
            if (!allowed.count(entry.first)) throw std::runtime_error("Unknown option " + entry.first);
        if (mode == "extract") {
            extract(args);
        } else if (mode == "fold") {
            fold(args);
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Zebra ch0 clock ERROR: " << error.what() << '\n';
        return 1;
    }
}
