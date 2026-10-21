#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

int main(int argc, char** argv) {
    try {
        if (argc != 4)
            throw std::runtime_error("Usage: extract_native_trigger_timestamps DECODED OUTPUT_DIRECTORY EXPECTED_WINDOWS");
        const fs::path inputPath(argv[1]), outputPath(argv[2]);
        const auto expected = std::stoull(argv[3]);
        if (fs::exists(outputPath)) throw std::runtime_error("Refusing existing output directory");
        const auto sizeBefore = fs::file_size(inputPath);
        const auto modifiedBefore = fs::last_write_time(inputPath);
        std::ifstream input(inputPath);
        if (!input) throw std::runtime_error("Cannot open decoded input");
        fs::create_directories(outputPath);
        std::ofstream binary(outputPath / "trigger_ticks_5ns.bin", std::ios::binary);
        std::ofstream samples(outputPath / "timestamp_samples.tsv");
        if (!binary || !samples) throw std::runtime_error("Cannot open output files");
        samples << "window_index\tdecoded_line\toriginal_seconds\tcoarse_tick_5ns\n";
        bool selected = false, anchorSeen = false;
        std::uint64_t count = 0, headers = 0, board0 = 0, lineNumber = 0;
        std::uint64_t first = 0, previous = 0;
        long double maxLatticeError = 0;
        std::string line;
        while (std::getline(input, line)) {
            ++lineNumber;
            if (line.rfind("TDC ", 0) == 0) {
                if (selected && !anchorSeen) throw std::runtime_error("Selected header without trigger timestamp");
                if (line.rfind("TDC 0 of ", 0) == 0) ++board0;
                selected = line.rfind("TDC 1 of ", 0) == 0;
                anchorSeen = false;
                if (selected) {
                    const auto location = line.rfind(" TDC ");
                    if (location == std::string::npos) throw std::runtime_error("Missing hardware identity");
                    std::istringstream label(line.substr(location + 5));
                    std::string hardware;
                    label >> hardware;
                    if (hardware != "TDC_001049")
                        throw std::runtime_error("Selected ordinal is not TDC_001049");
                    ++headers;
                }
                continue;
            }
            if (!selected || line.rfind("GlobalTriggerTime ", 0) != 0) continue;
            if (anchorSeen) throw std::runtime_error("Duplicate trigger timestamp in selected header");
            anchorSeen = true;
            const auto stamp = line.substr(18);
            const auto decimal = stamp.find('.');
            if (decimal == std::string::npos || stamp.size() - decimal - 1 < 9)
                throw std::runtime_error("At least nine timestamp decimals required");
            const long double ns = std::stold(stamp) * 1e9L;
            if (!std::isfinite(ns) || ns < 0) throw std::runtime_error("Invalid trigger timestamp");
            const auto tick = static_cast<std::uint64_t>(std::llround(ns / 5));
            const auto error = std::abs(ns - static_cast<long double>(tick) * 5);
            maxLatticeError = std::max(maxLatticeError, error);
            if (error > .65L) throw std::runtime_error("Timestamp does not match native 5 ns lattice");
            if (count && tick <= previous) throw std::runtime_error("Duplicate/backwards trigger timestamp");
            if (!count) first = tick;
            previous = tick;
            ++count;
            char bytes[8];
            for (int index = 0; index < 8; ++index)
                bytes[index] = static_cast<char>((tick >> (8 * index)) & 0xff);
            binary.write(bytes, 8);
            if (count == 1 || count % 1000 == 0 || count == expected)
                samples << count << '\t' << lineNumber << '\t' << stamp << '\t' << tick << '\n';
        }
        if (!input.eof()) throw std::runtime_error("Decoded input read error");
        if (selected && !anchorSeen) throw std::runtime_error("Final selected header has no timestamp");
        if (count != headers || count != expected || count < 2)
            throw std::runtime_error("Trigger count differs from independently validated selected-board count");
        binary.close();
        samples.close();
        if (!binary || !samples) throw std::runtime_error("Output write failure");
        if (sizeBefore != fs::file_size(inputPath) || modifiedBefore != fs::last_write_time(inputPath))
            throw std::runtime_error("Decoded input changed during extraction");
        std::ofstream summary(outputPath / "extraction_summary.txt");
        summary << std::setprecision(15)
                << "status=PASS\nboard=TDC_001049\ntdc_ordinal=1\ncoarse_quantum_ns=5\n"
                << "binary_format=little_endian_uint64_coarse_ticks\nselection=all_selected_board_headers_no_edge_or_reference_cuts\n"
                << "selected_board_headers=" << headers << "\ntrigger_timestamps=" << count
                << "\nboard0_headers=" << board0 << "\ndecoded_lines=" << lineNumber
                << "\nfirst_tick=" << first << "\nlast_tick=" << previous
                << "\nfirst_seconds=" << static_cast<long double>(first) * 5e-9L
                << "\nlast_seconds=" << static_cast<long double>(previous) * 5e-9L
                << "\nspan_seconds=" << static_cast<long double>(previous - first) * 5e-9L
                << "\nmax_lattice_rounding_error_ns=" << maxLatticeError
                << "\ninput_bytes=" << sizeBefore << '\n';
        summary.close();
        if (!summary) throw std::runtime_error("Summary write failure");
        std::cout << "PASS: extracted all " << count << " native board1049 trigger timestamps\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
