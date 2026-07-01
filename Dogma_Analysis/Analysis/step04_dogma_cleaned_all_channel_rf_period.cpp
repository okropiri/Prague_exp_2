#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <fcntl.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unistd.h>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr int kChannelCount = 32;
constexpr double kTwoPi = 6.28318530717958647692;
constexpr std::size_t kCycleResidualStatsLimit = 500000;
constexpr std::size_t kCycleResidualPointLimit = 500000;

struct Config {
    std::string inputPath;
    std::string outputRoot;
    std::string runKey;
    std::string ch0RefDirName = "Ch0_ref_Rates";
    std::string foldedDirName = "Folded_RF";
    std::string folded3xDirName = "Folded_RF_3x";
    std::string scanDirName = "RF_period_scan";
    double ch0RefMinNs = -6000.0;
    double ch0RefMaxNs = 6000.0;
    int ch0RefBins = 12000;
    double triggerRefMinNs = -6000.0;
    double triggerRefMaxNs = 6000.0;
    int triggerRefBins = 2400;
    double totMinNs = 0.0;
    double totMaxNs = 128.0;
    int totBins = 1280;
    int scoreChannel = 2;
    double scoreTotMinNs = 0.0;
    double scoreTotMaxNs = 128.0;
    double initialPeriodMinNs = 37.0;
    double initialPeriodMaxNs = 41.0;
    double initialStepNs = 0.05;
    int refineRounds = 3;
    int refineHalfSpanSteps = 25;
    double refineFactor = 10.0;
    double phaseBinWidthNs = 0.25;
    double peakWindowNs = 8.0;
    std::uint64_t minSelectedPulses = 1000;
    std::uint64_t scoreStride = 1;
    double ch0ValidRiseMinNs = -410.0;
    double ch0ValidRiseMaxNs = -395.0;
    double ch0ValidTotMinNs = 16.5;
    double ch0ValidTotMaxNs = 19.5;
    std::uint64_t segmentBlockScorePulses = 20000;
    int segmentPersistentBlocks = 2;
    double segmentPhaseJumpNs = 3.0;
    double segmentGapFactor = 100.0;
    double segmentMinGapSeconds = 1.0;
    bool singleRfSegment = false;
    bool skipCycleResidualDiagnostics = false;
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

struct StoredPulse {
    float timeNs = 0.0F;
    float totNs = 0.0F;
};

struct WindowState {
    bool hasValue = false;
    std::uint64_t windowIndex = 0;
    bool hasGlobalTriggerSeconds = false;
    double globalTriggerSeconds = 0.0;
    std::array<std::vector<PulseHit>, kChannelCount> hitsByChannel;
};

struct ScoredPulse {
    float timeNs = 0.0F;
    double globalTimeSeconds = 0.0;
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

struct PhaseBlock {
    std::size_t startPulseIndex = 0;
    std::size_t endPulseIndex = 0;
    double startGlobalTimeSeconds = 0.0;
    double endGlobalTimeSeconds = 0.0;
    CandidateMetrics fit;
};

struct PhaseSegment {
    int segmentId = 0;
    std::size_t startPulseIndex = 0;
    std::size_t endPulseIndex = 0;
    double startGlobalTimeSeconds = 0.0;
    double endGlobalTimeSeconds = 0.0;
    double phaseOriginNs = 0.0;
    std::uint64_t scorePulses = 0;
    double phaseShiftFromPreviousNs = 0.0;
    std::string reason = "initial";
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
    std::uint64_t scorePulsesBeforeStride = 0;
    std::uint64_t scorePulsesStored = 0;
    std::uint64_t ch0RefAcceptedPulses = 0;
    std::uint64_t foldedAcceptedPulses = 0;
    std::uint64_t folded3xAcceptedPulses = 0;
    std::uint64_t ch0RefRejectedOutsideRange = 0;
    std::array<std::uint64_t, kChannelCount> pulseCountByChannel{};
    std::array<std::uint64_t, kChannelCount> ch0RefAcceptedByChannel{};
    std::array<std::uint64_t, kChannelCount> foldedAcceptedByChannel{};
    std::array<std::uint64_t, kChannelCount> folded3xAcceptedByChannel{};
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

void drop_file_cache_hint(const std::string& path) {
#ifdef POSIX_FADV_DONTNEED
    const int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) {
        return;
    }
    (void)::posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED);
    ::close(fd);
#else
    (void)path;
#endif
}

double positive_mod(double value, double period) {
    const double wrapped = std::fmod(value, period);
    return wrapped < 0.0 ? wrapped + period : wrapped;
}

double wrap_centered(double value, double period) {
    return positive_mod(value + 0.5 * period, period) - 0.5 * period;
}

struct PhaseSeed {
    double centerNs = 0.0;
    std::uint64_t windowCount = 0;
};

std::vector<PhaseSeed> build_phase_seeds(const Config& config,
                                         const std::vector<std::uint64_t>& phaseCounts,
                                         double periodNs,
                                         double phaseBinWidthNs) {
    if (phaseCounts.empty()) {
        return {};
    }

    const int phaseBins = static_cast<int>(phaseCounts.size());
    std::vector<std::uint64_t> windowCounts(static_cast<std::size_t>(phaseBins), 0);
    for (int centerIndex = 0; centerIndex < phaseBins; ++centerIndex) {
        std::uint64_t windowCount = 0;
        for (int phaseIndex = 0; phaseIndex < phaseBins; ++phaseIndex) {
            const double centerNs = (static_cast<double>(centerIndex) + 0.5) * phaseBinWidthNs;
            const double phaseNs = (static_cast<double>(phaseIndex) + 0.5) * phaseBinWidthNs;
            const double distanceNs = std::abs(wrap_centered(phaseNs - centerNs, periodNs));
            if (distanceNs <= config.peakWindowNs) {
                windowCount += phaseCounts[static_cast<std::size_t>(phaseIndex)];
            }
        }
        windowCounts[static_cast<std::size_t>(centerIndex)] = windowCount;
    }

    std::vector<int> rankedIndices(static_cast<std::size_t>(phaseBins), 0);
    for (int phaseIndex = 0; phaseIndex < phaseBins; ++phaseIndex) {
        rankedIndices[static_cast<std::size_t>(phaseIndex)] = phaseIndex;
    }
    std::stable_sort(
        rankedIndices.begin(),
        rankedIndices.end(),
        [&](int leftIndex, int rightIndex) {
            const std::uint64_t leftWindowCount = windowCounts[static_cast<std::size_t>(leftIndex)];
            const std::uint64_t rightWindowCount = windowCounts[static_cast<std::size_t>(rightIndex)];
            if (leftWindowCount != rightWindowCount) {
                return leftWindowCount > rightWindowCount;
            }
            return phaseCounts[static_cast<std::size_t>(leftIndex)] > phaseCounts[static_cast<std::size_t>(rightIndex)];
        }
    );

    std::vector<PhaseSeed> seeds;
    seeds.reserve(4);
    for (int phaseIndex : rankedIndices) {
        const double centerNs = (static_cast<double>(phaseIndex) + 0.5) * phaseBinWidthNs;
        bool overlapsExistingSeed = false;
        for (const PhaseSeed& seed : seeds) {
            const double distanceNs = std::abs(wrap_centered(centerNs - seed.centerNs, periodNs));
            if (distanceNs < config.peakWindowNs) {
                overlapsExistingSeed = true;
                break;
            }
        }
        if (overlapsExistingSeed) {
            continue;
        }
        seeds.push_back(PhaseSeed{centerNs, windowCounts[static_cast<std::size_t>(phaseIndex)]});
        if (seeds.size() >= 4) {
            break;
        }
    }

    if (seeds.empty()) {
        const int bestIndex = rankedIndices.empty() ? 0 : rankedIndices.front();
        seeds.push_back(PhaseSeed{
            (static_cast<double>(bestIndex) + 0.5) * phaseBinWidthNs,
            windowCounts[static_cast<std::size_t>(bestIndex)]
        });
    }
    return seeds;
}

CandidateMetrics evaluate_seeded_candidate_range(const Config& config,
                                                 const std::vector<ScoredPulse>& scorePulses,
                                                 std::size_t beginIndex,
                                                 std::size_t endIndex,
                                                 double periodNs,
                                                 int roundIndex,
                                                 double seedCenterNs,
                                                 std::uint64_t seedWindowCount,
                                                 double coherence) {
    CandidateMetrics metrics;
    metrics.periodNs = periodNs;
    metrics.roundIndex = roundIndex;
    metrics.peakCenterNs = seedCenterNs;
    metrics.peakHeight = seedWindowCount;
    metrics.coherence = coherence;
    if (beginIndex >= endIndex || endIndex > scorePulses.size()) {
        return metrics;
    }

    double phaseOriginNs = seedCenterNs;
    for (int iteration = 0; iteration < 2; ++iteration) {
        double residualSum = 0.0;
        std::uint64_t selected = 0;
        for (std::size_t pulseIndex = beginIndex; pulseIndex < endIndex; ++pulseIndex) {
            const ScoredPulse& pulse = scorePulses[pulseIndex];
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
    for (std::size_t pulseIndex = beginIndex; pulseIndex < endIndex; ++pulseIndex) {
        const ScoredPulse& pulse = scorePulses[pulseIndex];
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
    metrics.selectedFraction = selectedDouble / static_cast<double>(endIndex - beginIndex);
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

CandidateMetrics evaluate_seeded_candidate(const Config& config,
                                           const std::vector<ScoredPulse>& scorePulses,
                                           double periodNs,
                                           int roundIndex,
                                           double seedCenterNs,
                                           std::uint64_t seedWindowCount,
                                           double coherence) {
    return evaluate_seeded_candidate_range(
        config,
        scorePulses,
        0,
        scorePulses.size(),
        periodNs,
        roundIndex,
        seedCenterNs,
        seedWindowCount,
        coherence
    );
}

bool is_better_candidate(const CandidateMetrics& candidate,
                         const CandidateMetrics& best) {
    if (candidate.valid != best.valid) {
        return candidate.valid;
    }
    if (candidate.valid && best.valid) {
        if (candidate.merit != best.merit) {
            return candidate.merit > best.merit;
        }
        if (candidate.selectedFraction != best.selectedFraction) {
            return candidate.selectedFraction > best.selectedFraction;
        }
        return candidate.sigmaNs < best.sigmaNs;
    }
    if (candidate.selectedPulses != best.selectedPulses) {
        return candidate.selectedPulses > best.selectedPulses;
    }
    if (candidate.peakHeight != best.peakHeight) {
        return candidate.peakHeight > best.peakHeight;
    }
    return candidate.peakCenterNs < best.peakCenterNs;
}

CandidateMetrics fit_phase_origin_for_range(const Config& config,
                                            const std::vector<ScoredPulse>& scorePulses,
                                            std::size_t beginIndex,
                                            std::size_t endIndex,
                                            double periodNs) {
    if (beginIndex >= endIndex || endIndex > scorePulses.size()) {
        return CandidateMetrics{};
    }

    const int phaseBins = std::max(16, static_cast<int>(std::llround(periodNs / config.phaseBinWidthNs)));
    std::vector<std::uint64_t> phaseCounts(static_cast<std::size_t>(phaseBins), 0);
    for (std::size_t pulseIndex = beginIndex; pulseIndex < endIndex; ++pulseIndex) {
        const double phaseNs = positive_mod(static_cast<double>(scorePulses[pulseIndex].timeNs), periodNs);
        const int phaseIndex = std::min(
            phaseBins - 1,
            static_cast<int>(std::floor(phaseNs / config.phaseBinWidthNs))
        );
        phaseCounts[static_cast<std::size_t>(phaseIndex)] += 1;
    }

    CandidateMetrics bestMetrics;
    const std::vector<PhaseSeed> seeds = build_phase_seeds(config, phaseCounts, periodNs, config.phaseBinWidthNs);
    for (const PhaseSeed& seed : seeds) {
        const CandidateMetrics candidate = evaluate_seeded_candidate_range(
            config,
            scorePulses,
            beginIndex,
            endIndex,
            periodNs,
            0,
            seed.centerNs,
            seed.windowCount,
            0.0
        );
        if (is_better_candidate(candidate, bestMetrics)) {
            bestMetrics = candidate;
        }
    }
    return bestMetrics;
}

double circular_mean_phase_origin(const std::vector<PhaseBlock>& blocks,
                                  std::size_t beginBlock,
                                  std::size_t endBlock,
                                  double periodNs,
                                  double fallbackPhaseOriginNs) {
    double cosineSum = 0.0;
    double sineSum = 0.0;
    const double angleScale = kTwoPi / periodNs;
    bool hasValue = false;
    for (std::size_t blockIndex = beginBlock; blockIndex < endBlock && blockIndex < blocks.size(); ++blockIndex) {
        const PhaseBlock& block = blocks[blockIndex];
        if (!block.fit.valid) {
            continue;
        }
        const double weight = std::max(1.0, static_cast<double>(block.fit.selectedPulses));
        const double angle = block.fit.phaseOriginNs * angleScale;
        cosineSum += weight * std::cos(angle);
        sineSum += weight * std::sin(angle);
        hasValue = true;
    }
    if (!hasValue || (cosineSum == 0.0 && sineSum == 0.0)) {
        return positive_mod(fallbackPhaseOriginNs, periodNs);
    }
    return positive_mod(std::atan2(sineSum, cosineSum) / angleScale, periodNs);
}

double compute_gap_threshold_seconds(const Config& config,
                                     const std::vector<ScoredPulse>& scorePulses) {
    if (scorePulses.size() < 2) {
        return config.segmentMinGapSeconds;
    }

    const std::size_t sampleLimit = 100000;
    const std::size_t stride = std::max<std::size_t>(1, (scorePulses.size() - 1) / sampleLimit);
    std::vector<double> gaps;
    gaps.reserve(std::min(sampleLimit, scorePulses.size() - 1));
    for (std::size_t pulseIndex = 1; pulseIndex < scorePulses.size(); pulseIndex += stride) {
        const double gapSeconds = scorePulses[pulseIndex].globalTimeSeconds - scorePulses[pulseIndex - 1].globalTimeSeconds;
        if (gapSeconds > 0.0 && std::isfinite(gapSeconds)) {
            gaps.push_back(gapSeconds);
        }
    }
    if (gaps.empty()) {
        return config.segmentMinGapSeconds;
    }

    const std::size_t medianIndex = gaps.size() / 2;
    std::nth_element(gaps.begin(), gaps.begin() + static_cast<std::ptrdiff_t>(medianIndex), gaps.end());
    const double medianGapSeconds = gaps[medianIndex];
    return std::max(config.segmentMinGapSeconds, config.segmentGapFactor * medianGapSeconds);
}

std::vector<PhaseBlock> build_phase_blocks(const Config& config,
                                           const std::vector<ScoredPulse>& scorePulses,
                                           double periodNs) {
    std::vector<PhaseBlock> blocks;
    if (scorePulses.empty()) {
        return blocks;
    }

    const double gapThresholdSeconds = compute_gap_threshold_seconds(config, scorePulses);
    std::size_t startIndex = 0;
    while (startIndex < scorePulses.size()) {
        std::size_t endIndex = std::min(scorePulses.size(), startIndex + config.segmentBlockScorePulses);
        for (std::size_t pulseIndex = startIndex + 1; pulseIndex < endIndex; ++pulseIndex) {
            const double gapSeconds = scorePulses[pulseIndex].globalTimeSeconds - scorePulses[pulseIndex - 1].globalTimeSeconds;
            if (gapSeconds > gapThresholdSeconds) {
                endIndex = pulseIndex;
                break;
            }
        }
        if (endIndex <= startIndex) {
            endIndex = std::min(scorePulses.size(), startIndex + 1);
        }

        PhaseBlock block;
        block.startPulseIndex = startIndex;
        block.endPulseIndex = endIndex;
        block.startGlobalTimeSeconds = scorePulses[startIndex].globalTimeSeconds;
        block.endGlobalTimeSeconds = scorePulses[endIndex - 1].globalTimeSeconds;
        block.fit = fit_phase_origin_for_range(config, scorePulses, startIndex, endIndex, periodNs);
        blocks.push_back(block);
        startIndex = endIndex;
    }

    return blocks;
}

PhaseSegment fit_phase_segment(const Config& config,
                               const std::vector<ScoredPulse>& scorePulses,
                               const std::vector<PhaseBlock>& blocks,
                               std::size_t beginBlock,
                               std::size_t endBlock,
                               const CandidateMetrics& best,
                               int segmentId,
                               const std::string& reason,
                               double previousPhaseOriginNs) {
    PhaseSegment segment;
    segment.segmentId = segmentId;
    segment.reason = reason;
    segment.phaseOriginNs = best.phaseOriginNs;
    if (beginBlock >= endBlock || endBlock > blocks.size()) {
        return segment;
    }

    segment.startPulseIndex = blocks[beginBlock].startPulseIndex;
    segment.endPulseIndex = blocks[endBlock - 1].endPulseIndex;
    segment.scorePulses = static_cast<std::uint64_t>(segment.endPulseIndex - segment.startPulseIndex);
    if (segment.startPulseIndex < scorePulses.size() && segment.endPulseIndex > segment.startPulseIndex) {
        segment.startGlobalTimeSeconds = scorePulses[segment.startPulseIndex].globalTimeSeconds;
        segment.endGlobalTimeSeconds = scorePulses[segment.endPulseIndex - 1].globalTimeSeconds;
    }

    const CandidateMetrics segmentFit = fit_phase_origin_for_range(
        config,
        scorePulses,
        segment.startPulseIndex,
        segment.endPulseIndex,
        best.periodNs
    );
    if (segmentFit.valid) {
        segment.phaseOriginNs = segmentFit.phaseOriginNs;
    }
    if (segmentId > 0) {
        segment.phaseShiftFromPreviousNs = std::abs(
            wrap_centered(segment.phaseOriginNs - previousPhaseOriginNs, best.periodNs)
        );
    }
    return segment;
}

std::vector<PhaseSegment> build_phase_segments(const Config& config,
                                               const std::vector<ScoredPulse>& scorePulses,
                                               const std::vector<PhaseBlock>& blocks,
                                               const CandidateMetrics& best) {
    std::vector<PhaseSegment> segments;
    if (scorePulses.empty()) {
        return segments;
    }
    if (blocks.empty()) {
        PhaseSegment segment;
        segment.segmentId = 0;
        segment.startPulseIndex = 0;
        segment.endPulseIndex = scorePulses.size();
        segment.startGlobalTimeSeconds = scorePulses.front().globalTimeSeconds;
        segment.endGlobalTimeSeconds = scorePulses.back().globalTimeSeconds;
        segment.phaseOriginNs = best.phaseOriginNs;
        segment.scorePulses = static_cast<std::uint64_t>(scorePulses.size());
        segments.push_back(segment);
        return segments;
    }

    std::size_t segmentStartBlock = 0;
    std::string currentReason = "initial";
    std::size_t blockIndex = 1;
    while (blockIndex < blocks.size()) {
        const PhaseBlock& currentBlock = blocks[blockIndex];
        if (!currentBlock.fit.valid) {
            ++blockIndex;
            continue;
        }

        const double baselinePhaseOriginNs = circular_mean_phase_origin(
            blocks,
            segmentStartBlock,
            blockIndex,
            best.periodNs,
            best.phaseOriginNs
        );
        const double phaseShiftNs = std::abs(
            wrap_centered(currentBlock.fit.phaseOriginNs - baselinePhaseOriginNs, best.periodNs)
        );
        if (phaseShiftNs < config.segmentPhaseJumpNs) {
            ++blockIndex;
            continue;
        }

        int confirmedBlocks = 0;
        std::size_t lookahead = blockIndex;
        double candidatePhaseOriginNs = currentBlock.fit.phaseOriginNs;
        while (lookahead < blocks.size() && confirmedBlocks < config.segmentPersistentBlocks) {
            const PhaseBlock& candidateBlock = blocks[lookahead];
            if (!candidateBlock.fit.valid) {
                break;
            }
            const double shiftFromBaselineNs = std::abs(
                wrap_centered(candidateBlock.fit.phaseOriginNs - baselinePhaseOriginNs, best.periodNs)
            );
            const double shiftFromCandidateNs = std::abs(
                wrap_centered(candidateBlock.fit.phaseOriginNs - candidatePhaseOriginNs, best.periodNs)
            );
            if (shiftFromBaselineNs < config.segmentPhaseJumpNs) {
                break;
            }
            if (confirmedBlocks > 0 && shiftFromCandidateNs > config.segmentPhaseJumpNs) {
                break;
            }
            candidatePhaseOriginNs = circular_mean_phase_origin(
                blocks,
                blockIndex,
                lookahead + 1,
                best.periodNs,
                candidatePhaseOriginNs
            );
            ++confirmedBlocks;
            ++lookahead;
        }

        const bool strongTerminalShift = confirmedBlocks >= 1
            && lookahead >= blocks.size()
            && phaseShiftNs >= (2.0 * config.segmentPhaseJumpNs);
        if (confirmedBlocks >= config.segmentPersistentBlocks || strongTerminalShift) {
            const double previousPhaseOriginNs = segments.empty() ? best.phaseOriginNs : segments.back().phaseOriginNs;
            segments.push_back(fit_phase_segment(
                config,
                scorePulses,
                blocks,
                segmentStartBlock,
                blockIndex,
                best,
                static_cast<int>(segments.size()),
                currentReason,
                previousPhaseOriginNs
            ));
            segmentStartBlock = blockIndex;
            currentReason = "persistent_phase_jump";
            continue;
        }

        ++blockIndex;
    }

    const double previousPhaseOriginNs = segments.empty() ? best.phaseOriginNs : segments.back().phaseOriginNs;
    segments.push_back(fit_phase_segment(
        config,
        scorePulses,
        blocks,
        segmentStartBlock,
        blocks.size(),
        best,
        static_cast<int>(segments.size()),
        currentReason,
        previousPhaseOriginNs
    ));
    return segments;
}

PhaseSegment build_single_phase_segment(const std::vector<ScoredPulse>& scorePulses,
                                        const CandidateMetrics& best) {
    PhaseSegment segment;
    segment.segmentId = 0;
    segment.startPulseIndex = 0;
    segment.endPulseIndex = scorePulses.size();
    segment.phaseOriginNs = best.phaseOriginNs;
    segment.scorePulses = static_cast<std::uint64_t>(scorePulses.size());
    segment.reason = "single_segment";
    if (!scorePulses.empty()) {
        segment.startGlobalTimeSeconds = scorePulses.front().globalTimeSeconds;
        segment.endGlobalTimeSeconds = scorePulses.back().globalTimeSeconds;
    }
    return segment;
}

const PhaseSegment* find_phase_segment_for_time(const std::vector<PhaseSegment>& segments,
                                                double globalTimeSeconds) {
    if (segments.empty() || !std::isfinite(globalTimeSeconds)) {
        return nullptr;
    }

    const auto iterator = std::upper_bound(
        segments.begin(),
        segments.end(),
        globalTimeSeconds,
        [](double timeSeconds, const PhaseSegment& segment) {
            return timeSeconds < segment.startGlobalTimeSeconds;
        }
    );
    if (iterator == segments.begin()) {
        return &segments.front();
    }
    return &*(iterator - 1);
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

void write_histogram_1d(const fs::path& path, const Histogram1D& histogram) {
    std::ofstream output(path);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write histogram file: " + path.string());
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
        const std::uint64_t count = histogram.storage[index];
        if (count == 0) {
            continue;
        }
        output << index << ' ' << count << '\n';
    }
}

void write_histogram_2d(const fs::path& path, const Histogram2D& histogram) {
    std::ofstream output(path);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write histogram file: " + path.string());
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
            const std::uint64_t count = histogram.storage[flatIndex];
            if (count == 0) {
                continue;
            }
            output << xIndex << ' ' << yIndex << ' ' << count << '\n';
        }
    }
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
        } else if (arg == "--ch0ref-dir-name" && index + 1 < argc) {
            config.ch0RefDirName = argv[++index];
        } else if (arg == "--folded-dir-name" && index + 1 < argc) {
            config.foldedDirName = argv[++index];
        } else if (arg == "--folded3x-dir-name" && index + 1 < argc) {
            config.folded3xDirName = argv[++index];
        } else if (arg == "--scan-dir-name" && index + 1 < argc) {
            config.scanDirName = argv[++index];
        } else if (arg == "--ch0ref-min-ns" && index + 1 < argc) {
            config.ch0RefMinNs = std::stod(argv[++index]);
        } else if (arg == "--ch0ref-max-ns" && index + 1 < argc) {
            config.ch0RefMaxNs = std::stod(argv[++index]);
        } else if (arg == "--ch0ref-bins" && index + 1 < argc) {
            config.ch0RefBins = std::stoi(argv[++index]);
        } else if (arg == "--trigger-ref-min-ns" && index + 1 < argc) {
            config.triggerRefMinNs = std::stod(argv[++index]);
        } else if (arg == "--trigger-ref-max-ns" && index + 1 < argc) {
            config.triggerRefMaxNs = std::stod(argv[++index]);
        } else if (arg == "--trigger-ref-bins" && index + 1 < argc) {
            config.triggerRefBins = std::stoi(argv[++index]);
        } else if (arg == "--tot-min-ns" && index + 1 < argc) {
            config.totMinNs = std::stod(argv[++index]);
        } else if (arg == "--tot-max-ns" && index + 1 < argc) {
            config.totMaxNs = std::stod(argv[++index]);
        } else if (arg == "--tot-bins" && index + 1 < argc) {
            config.totBins = std::stoi(argv[++index]);
        } else if (arg == "--score-channel" && index + 1 < argc) {
            config.scoreChannel = std::stoi(argv[++index]);
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
        } else if (arg == "--ch0-valid-rise-min-ns" && index + 1 < argc) {
            config.ch0ValidRiseMinNs = std::stod(argv[++index]);
        } else if (arg == "--ch0-valid-rise-max-ns" && index + 1 < argc) {
            config.ch0ValidRiseMaxNs = std::stod(argv[++index]);
        } else if (arg == "--ch0-valid-tot-min-ns" && index + 1 < argc) {
            config.ch0ValidTotMinNs = std::stod(argv[++index]);
        } else if (arg == "--ch0-valid-tot-max-ns" && index + 1 < argc) {
            config.ch0ValidTotMaxNs = std::stod(argv[++index]);
        } else if (arg == "--segment-block-score-pulses" && index + 1 < argc) {
            config.segmentBlockScorePulses = static_cast<std::uint64_t>(std::stoull(argv[++index]));
        } else if (arg == "--segment-persistent-blocks" && index + 1 < argc) {
            config.segmentPersistentBlocks = std::stoi(argv[++index]);
        } else if (arg == "--segment-phase-jump-ns" && index + 1 < argc) {
            config.segmentPhaseJumpNs = std::stod(argv[++index]);
        } else if (arg == "--segment-gap-factor" && index + 1 < argc) {
            config.segmentGapFactor = std::stod(argv[++index]);
        } else if (arg == "--segment-min-gap-seconds" && index + 1 < argc) {
            config.segmentMinGapSeconds = std::stod(argv[++index]);
        } else if (arg == "--single-rf-segment") {
            config.singleRfSegment = true;
        } else if (arg == "--skip-cycle-residual-diagnostics") {
            config.skipCycleResidualDiagnostics = true;
        } else if (arg == "--help") {
            std::cout
                << "Usage: dogma_cleaned_all_channel_rf_period --input <cleaned_pulses.tsv> --output-root <Results_root>\n"
                << "  [--run-key <name>] [--ch0ref-dir-name Ch0_ref_Rates] [--folded-dir-name Folded_RF] [--folded3x-dir-name Folded_RF_3x] [--scan-dir-name RF_period_scan]\n"
                << "  [--ch0ref-min-ns -6000] [--ch0ref-max-ns 6000] [--ch0ref-bins 12000]\n"
                << "  [--trigger-ref-min-ns -6000] [--trigger-ref-max-ns 6000] [--trigger-ref-bins 2400]\n"
                << "  [--tot-min-ns 0] [--tot-max-ns 128] [--tot-bins 1280]\n"
                << "  [--score-channel 2] [--score-tot-min-ns 0] [--score-tot-max-ns 128]\n"
                << "  [--initial-period-min-ns 37.0] [--initial-period-max-ns 41.0] [--initial-step-ns 0.05]\n"
                << "  [--refine-rounds 3] [--refine-half-span-steps 5] [--refine-factor 10]\n"
                << "  [--phase-bin-width-ns 0.25] [--peak-window-ns 8] [--min-selected-pulses 1000] [--score-stride 1]\n"
                << "  [--ch0-valid-rise-min-ns -410] [--ch0-valid-rise-max-ns -395]\n"
                << "  [--ch0-valid-tot-min-ns 16.5] [--ch0-valid-tot-max-ns 19.5]\n"
                << "  [--segment-block-score-pulses 20000] [--segment-persistent-blocks 2]\n"
                << "  [--segment-phase-jump-ns 3.0] [--segment-gap-factor 100] [--segment-min-gap-seconds 1.0]\n"
                << "  [--single-rf-segment] [--skip-cycle-residual-diagnostics]\n";
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
    if (!(config.ch0RefMinNs < config.ch0RefMaxNs)) {
        throw std::runtime_error("ch0ref-min-ns must be smaller than ch0ref-max-ns");
    }
    if (config.ch0RefBins <= 0) {
        throw std::runtime_error("ch0ref-bins must be positive");
    }
    if (!(config.triggerRefMinNs < config.triggerRefMaxNs)) {
        throw std::runtime_error("trigger-ref-min-ns must be smaller than trigger-ref-max-ns");
    }
    if (config.triggerRefBins <= 0) {
        throw std::runtime_error("trigger-ref-bins must be positive");
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
    if (config.refineRounds <= 0 || config.refineHalfSpanSteps <= 0) {
        throw std::runtime_error("refinement parameters must be positive");
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
    if (config.minSelectedPulses == 0 || config.scoreStride == 0) {
        throw std::runtime_error("min-selected-pulses and score-stride must be positive");
    }
    if (!(config.ch0ValidRiseMinNs < config.ch0ValidRiseMaxNs)) {
        throw std::runtime_error("ch0-valid-rise-min-ns must be smaller than ch0-valid-rise-max-ns");
    }
    if (!(config.ch0ValidTotMinNs < config.ch0ValidTotMaxNs)) {
        throw std::runtime_error("ch0-valid-tot-min-ns must be smaller than ch0-valid-tot-max-ns");
    }
    if (config.segmentBlockScorePulses == 0) {
        throw std::runtime_error("segment-block-score-pulses must be positive");
    }
    if (config.segmentPersistentBlocks <= 0) {
        throw std::runtime_error("segment-persistent-blocks must be positive");
    }
    if (config.segmentPhaseJumpNs <= 0.0) {
        throw std::runtime_error("segment-phase-jump-ns must be positive");
    }
    if (config.segmentGapFactor <= 0.0) {
        throw std::runtime_error("segment-gap-factor must be positive");
    }
    if (config.segmentMinGapSeconds < 0.0) {
        throw std::runtime_error("segment-min-gap-seconds must be non-negative");
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

bool passes_score_tot_gate(const Config& config, double totNs) {
    return config.scoreTotMinNs <= totNs && totNs <= config.scoreTotMaxNs;
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

using SparseHistogram = std::array<std::vector<std::uint32_t>, kChannelCount>;

std::uint64_t sparse_key(int xBinZeroBased, int yBinZeroBased, int yBins) {
    return static_cast<std::uint64_t>(xBinZeroBased) * static_cast<std::uint64_t>(yBins)
         + static_cast<std::uint64_t>(yBinZeroBased);
}

void ensure_sparse_histogram_sizes(SparseHistogram& histograms, int xBins, int yBins, int storedChannel) {
    const std::size_t size = static_cast<std::size_t>(xBins) * static_cast<std::size_t>(yBins);
    for (int channel = 0; channel < kChannelCount; ++channel) {
        auto& histogram = histograms[static_cast<std::size_t>(channel)];
        if (channel == storedChannel) {
            histogram.assign(size, 0U);
        } else {
            histogram.clear();
        }
    }
}

void fill_sparse_histogram(SparseHistogram& histogram,
                           int channel,
                           double xValue,
                           int xBins,
                           double xMin,
                           double xMax,
                           double yValue,
                           int yBins,
                           double yMin,
                           double yMax) {
    auto& storage = histogram[static_cast<std::size_t>(channel)];
    if (storage.empty()) {
        return;
    }
    const int xBin = histogram_bin_index(xValue, xBins, xMin, xMax);
    const int yBin = histogram_bin_index(yValue, yBins, yMin, yMax);
    if (xBin <= 0 || xBin > xBins || yBin <= 0 || yBin > yBins) {
        return;
    }
    const std::size_t key = static_cast<std::size_t>(sparse_key(xBin - 1, yBin - 1, yBins));
    storage[key] += 1U;
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

    for (const ScoredPulse& pulse : scorePulses) {
        const double timeNs = static_cast<double>(pulse.timeNs);
        const double phaseNs = positive_mod(timeNs, periodNs);
        const int phaseIndex = std::min(phaseBins - 1, static_cast<int>(std::floor(phaseNs / phaseBinWidthNs)));
        phaseCounts[static_cast<std::size_t>(phaseIndex)] += 1;
    }

    double sumCos = 0.0;
    double sumSin = 0.0;
    for (int phaseIndex = 0; phaseIndex < phaseBins; ++phaseIndex) {
        const double binCenterNs = (static_cast<double>(phaseIndex) + 0.5) * phaseBinWidthNs;
        const double angle = kTwoPi * binCenterNs / periodNs;
        const double count = static_cast<double>(phaseCounts[static_cast<std::size_t>(phaseIndex)]);
        sumCos += count * std::cos(angle);
        sumSin += count * std::sin(angle);
    }
    const double coherence = std::sqrt(sumCos * sumCos + sumSin * sumSin) / static_cast<double>(scorePulses.size());
    const std::vector<PhaseSeed> phaseSeeds = build_phase_seeds(config, phaseCounts, periodNs, phaseBinWidthNs);
    bool hasBestSeed = false;
    for (const PhaseSeed& seed : phaseSeeds) {
        const CandidateMetrics seededMetrics = evaluate_seeded_candidate(
            config,
            scorePulses,
            periodNs,
            roundIndex,
            seed.centerNs,
            seed.windowCount,
            coherence
        );
        if (!hasBestSeed || is_better_candidate(seededMetrics, metrics)) {
            metrics = seededMetrics;
            hasBestSeed = true;
        }
    }
    return metrics;
}

std::vector<CandidateMetrics> scan_candidates(const Config& config,
                                              const std::vector<ScoredPulse>& scorePulses) {
    std::vector<CandidateMetrics> rows;
    double periodMin = config.initialPeriodMinNs;
    double periodMax = config.initialPeriodMaxNs;
    double stepNs = config.initialStepNs;
    CandidateMetrics best;
    bool hasBest = false;

    for (int roundIndex = 0; roundIndex < config.refineRounds; ++roundIndex) {
        for (double periodNs = periodMin; periodNs <= periodMax + 0.5 * stepNs; periodNs += stepNs) {
            const CandidateMetrics candidate = evaluate_candidate(config, scorePulses, periodNs, roundIndex);
            rows.push_back(candidate);
            if (candidate.valid && (!hasBest || candidate.merit > best.merit)) {
                best = candidate;
                hasBest = true;
            }
        }
        if (!hasBest) {
            continue;
        }
        stepNs /= config.refineFactor;
        periodMin = best.periodNs - static_cast<double>(config.refineHalfSpanSteps) * stepNs;
        periodMax = best.periodNs + static_cast<double>(config.refineHalfSpanSteps) * stepNs;
    }
    return rows;
}

CandidateMetrics best_valid_candidate(const std::vector<CandidateMetrics>& rows) {
    CandidateMetrics best;
    bool hasBest = false;
    for (const CandidateMetrics& row : rows) {
        if (!row.valid) {
            continue;
        }
        if (!hasBest || row.merit > best.merit) {
            best = row;
            hasBest = true;
        }
    }
    if (!hasBest) {
        throw std::runtime_error("RF period scan did not find a valid candidate");
    }
    return best;
}

Histogram1D build_best_phase_profile(const Config& config,
                                     const std::vector<ScoredPulse>& scorePulses,
                                     const CandidateMetrics& best) {
    const int phaseBins = std::max(16, static_cast<int>(std::llround(best.periodNs / config.phaseBinWidthNs)));
    Histogram1D histogram = make_histogram_1d(
        "cleaned_RfPhaseProfile",
        "RF-phase profile from cleaned ch0-referenced score-channel times",
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

Histogram2D make_best_phase_tot_histogram(const Config& config,
                                          const CandidateMetrics& best) {
    const int phaseBins = std::max(16, static_cast<int>(std::llround(best.periodNs / config.phaseBinWidthNs)));
    return make_histogram_2d(
        "cleaned_RfPhaseTot",
        "Score-channel ToT vs recovered RF phase from cleaned ch0-referenced times",
        phaseBins,
        0.0,
        best.periodNs,
        "RF phase (ns)",
        config.totBins,
        config.totMinNs,
        config.totMaxNs,
        "ToT (ns)"
    );
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
    // Long runs can touch millions of RF cycles. Keep this diagnostic bounded so
    // it cannot dominate the memory footprint after the RF solution is found.
    const std::size_t sourceStride = std::max<std::size_t>(
        1,
        (scorePulses.size() + kCycleResidualStatsLimit - 1) / kCycleResidualStatsLimit
    );
    for (std::size_t pulseIndex = 0; pulseIndex < scorePulses.size(); pulseIndex += sourceStride) {
        const ScoredPulse& pulse = scorePulses[pulseIndex];
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
    if (scorePulses.empty()) {
        return points;
    }

    // The diagnostic scatter does not need every selected pulse on the largest runs.
    const std::size_t sourceStride = std::max<std::size_t>(
        1,
        (scorePulses.size() + kCycleResidualPointLimit - 1) / kCycleResidualPointLimit
    );
    points.reserve(std::min(kCycleResidualPointLimit, scorePulses.size()));
    for (std::size_t pulseIndex = 0; pulseIndex < scorePulses.size(); pulseIndex += sourceStride) {
        const ScoredPulse& pulse = scorePulses[pulseIndex];
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

void write_scan_file(const fs::path& outputPath,
                     const Config& config,
                     const std::vector<CandidateMetrics>& rows) {
    std::ofstream output(outputPath);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write RF scan file: " + outputPath.string());
    }
    output << "# initial_period_range_ns=[" << format_double(config.initialPeriodMinNs) << ", "
           << format_double(config.initialPeriodMaxNs) << "]\n";
    output << "# initial_step_ns=" << format_double(config.initialStepNs) << '\n';
    output << "# refine_rounds=" << config.refineRounds << '\n';
    output << "# refine_half_span_steps=" << config.refineHalfSpanSteps << '\n';
    output << "# refine_factor=" << format_double(config.refineFactor) << '\n';
    output << "# phase_bin_width_ns=" << format_double(config.phaseBinWidthNs) << '\n';
    output << "# peak_window_ns=" << format_double(config.peakWindowNs) << '\n';
    output << "# score_tot_gate_ns=[" << format_double(config.scoreTotMinNs) << ", "
           << format_double(config.scoreTotMaxNs) << "]\n";
    output << "# columns: round_index period_ns phase_origin_ns peak_center_ns peak_height selected_pulses selected_fraction sigma_ns mean_residual_ns coherence drift_slope_ns_per_cycle drift_intercept_ns merit valid\n";
    output << "# peak_center_ns is the best seed-window center searched across the full folded phase; peak_height is the pulse count in that seed window before recentering.\n";
    for (const CandidateMetrics& row : rows) {
        output << row.roundIndex << ' '
               << format_double(row.periodNs) << ' '
               << format_double(row.phaseOriginNs) << ' '
               << format_double(row.peakCenterNs) << ' '
               << row.peakHeight << ' '
               << row.selectedPulses << ' '
               << format_double(row.selectedFraction) << ' '
               << format_double(row.sigmaNs) << ' '
               << format_double(row.meanResidualNs) << ' '
               << format_double(row.coherence) << ' '
               << format_double(row.driftSlopeNsPerCycle) << ' '
               << format_double(row.driftInterceptNs) << ' '
               << format_double(row.merit) << ' '
               << (row.valid ? "true" : "false") << '\n';
    }
}

void write_cycle_residual_file(const fs::path& outputPath,
                               const CandidateMetrics& best,
                               const std::vector<CycleResidualStats>& rows) {
    std::ofstream output(outputPath);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write cycle residual file: " + outputPath.string());
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

void write_cycle_residual_points_file(const fs::path& outputPath,
                                      const CandidateMetrics& best,
                                      const std::vector<CycleResidualPoint>& points) {
    std::ofstream output(outputPath);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write cycle residual points file: " + outputPath.string());
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

void finalize_window_for_scan(const WindowState& window,
                              const Config& config,
                              Stats& stats,
                              std::vector<ScoredPulse>& scorePulses,
                              std::uint64_t& scoreOrdinal) {
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

    const auto& scoreHits = window.hitsByChannel[static_cast<std::size_t>(config.scoreChannel)];
    for (const PulseHit& hit : scoreHits) {
        if (!passes_score_tot_gate(config, hit.totNs)) {
            continue;
        }
        if ((scoreOrdinal % config.scoreStride) == 0) {
            const double globalTimeSeconds = window.hasGlobalTriggerSeconds
                ? window.globalTriggerSeconds + hit.riseNs * 1.0e-9
                : hit.riseNs * 1.0e-9;
            scorePulses.push_back(ScoredPulse{static_cast<float>(hit.riseNs - referenceHit->riseNs), globalTimeSeconds});
            stats.scorePulsesStored += 1;
        }
        scoreOrdinal += 1;
        stats.scorePulsesBeforeStride += 1;
    }
}

using ChannelProfiles = std::array<std::vector<std::uint64_t>, kChannelCount>;
using ChannelPhaseTotHistograms = std::array<std::vector<std::uint64_t>, kChannelCount>;

void ensure_profile_sizes(ChannelProfiles& profiles, int bins) {
    for (auto& profile : profiles) {
        profile.assign(static_cast<std::size_t>(bins), 0);
    }
}

void ensure_phase_tot_sizes(ChannelPhaseTotHistograms& histograms, int xBins, int yBins, int storedChannel) {
    for (int channel = 0; channel < kChannelCount; ++channel) {
        auto& histogram = histograms[static_cast<std::size_t>(channel)];
        if (channel == storedChannel) {
            histogram.assign(static_cast<std::size_t>(xBins) * static_cast<std::size_t>(yBins), 0);
        } else {
            histogram.clear();
        }
    }
}

void fill_profile(ChannelProfiles& profiles,
                  std::vector<std::uint64_t>& totals,
                  int channel,
                  double value,
                  int bins,
                  double min,
                  double max) {
    const int bin = histogram_bin_index(value, bins, min, max);
    if (bin <= 0 || bin > bins) {
        return;
    }
    const std::size_t index = static_cast<std::size_t>(bin - 1);
    profiles[static_cast<std::size_t>(channel)][index] += 1;
    totals[index] += 1;
}

void fill_phase_tot(ChannelPhaseTotHistograms& histograms,
                    int channel,
                    double xValue,
                    int xBins,
                    double xMin,
                    double xMax,
                    double totValue,
                    int totBins,
                    double totMin,
                    double totMax) {
    auto& storage = histograms[static_cast<std::size_t>(channel)];
    if (storage.empty()) {
        return;
    }
    const int xBin = histogram_bin_index(xValue, xBins, xMin, xMax);
    const int yBin = histogram_bin_index(totValue, totBins, totMin, totMax);
    if (xBin <= 0 || xBin > xBins || yBin <= 0 || yBin > totBins) {
        return;
    }
    const std::size_t flatIndex = static_cast<std::size_t>(yBin - 1) * static_cast<std::size_t>(xBins)
                                + static_cast<std::size_t>(xBin - 1);
    storage[flatIndex] += 1;
}

void finalize_window_for_profiles(const WindowState& window,
                                  const Config& config,
                                  const CandidateMetrics& best,
                                  Histogram2D& bestPhaseTot,
                                  const std::vector<PhaseSegment>& phaseSegments,
                                  Stats& stats,
                                  ChannelProfiles& ch0RefProfiles,
                                  std::vector<std::uint64_t>& ch0RefTotals,
                                  ChannelProfiles& foldedProfiles,
                                  std::vector<std::uint64_t>& foldedTotals,
                                  ChannelPhaseTotHistograms& foldedPhaseTotHistograms,
                                  SparseHistogram& foldedPhaseCh0TimeHistograms,
                                  SparseHistogram& foldedPhaseTriggerTimeHistograms,
                                  ChannelProfiles& folded3xProfiles,
                                  std::vector<std::uint64_t>& folded3xTotals,
                                  ChannelPhaseTotHistograms& folded3xPhaseTotHistograms,
                                  SparseHistogram& folded3xPhaseCh0TimeHistograms,
                                  SparseHistogram& folded3xPhaseTriggerTimeHistograms) {
    if (!window.hasValue) {
        return;
    }

    const auto& ch0Hits = window.hitsByChannel[0];
    if (ch0Hits.empty()) {
        return;
    }

    const PulseHit* referenceHit = nullptr;
    for (const PulseHit& hit : ch0Hits) {
        if (hit.riseNs < config.ch0ValidRiseMinNs || hit.riseNs > config.ch0ValidRiseMaxNs) {
            continue;
        }
        if (hit.totNs < config.ch0ValidTotMinNs || hit.totNs > config.ch0ValidTotMaxNs) {
            continue;
        }
        referenceHit = &hit;
        break;
    }
    if (referenceHit == nullptr) {
        return;
    }

    for (const PulseHit& hit : window.hitsByChannel[static_cast<std::size_t>(config.scoreChannel)]) {
        const double scorePhaseNs = positive_mod((hit.riseNs - referenceHit->riseNs) - best.phaseOriginNs, best.periodNs);
        fill_histogram(bestPhaseTot, scorePhaseNs, hit.totNs);
    }

    const int foldedBins = static_cast<int>(foldedTotals.size());
    const int folded3xBins = static_cast<int>(folded3xTotals.size());
    const double period3xNs = 3.0 * best.periodNs;
    const PhaseSegment* phaseSegment = find_phase_segment_for_time(phaseSegments, window.globalTriggerSeconds);
    const double phaseOriginNs = phaseSegment == nullptr ? best.phaseOriginNs : phaseSegment->phaseOriginNs;

    auto fill_all_modes = [&](int channel, double relativeNs, double triggerRelativeNs, double totNs) {
        if (relativeNs < config.ch0RefMinNs || relativeNs >= config.ch0RefMaxNs) {
            stats.ch0RefRejectedOutsideRange += 1;
        } else {
            fill_profile(ch0RefProfiles, ch0RefTotals, channel, relativeNs, config.ch0RefBins, config.ch0RefMinNs, config.ch0RefMaxNs);
            stats.ch0RefAcceptedPulses += 1;
            stats.ch0RefAcceptedByChannel[static_cast<std::size_t>(channel)] += 1;
        }

        const double phaseNs = positive_mod(relativeNs - phaseOriginNs, best.periodNs);
        fill_profile(foldedProfiles, foldedTotals, channel, phaseNs, foldedBins, 0.0, best.periodNs);
        if (channel == config.scoreChannel) {
            fill_phase_tot(foldedPhaseTotHistograms, channel, phaseNs, foldedBins, 0.0, best.periodNs, totNs, config.totBins, config.totMinNs, config.totMaxNs);
            fill_sparse_histogram(
                foldedPhaseCh0TimeHistograms,
                channel,
                phaseNs,
                foldedBins,
                0.0,
                best.periodNs,
                relativeNs,
                config.ch0RefBins,
                config.ch0RefMinNs,
                config.ch0RefMaxNs
            );
            fill_sparse_histogram(
                foldedPhaseTriggerTimeHistograms,
                channel,
                phaseNs,
                foldedBins,
                0.0,
                best.periodNs,
                triggerRelativeNs,
                config.triggerRefBins,
                config.triggerRefMinNs,
                config.triggerRefMaxNs
            );
        }
        stats.foldedAcceptedPulses += 1;
        stats.foldedAcceptedByChannel[static_cast<std::size_t>(channel)] += 1;

        const double phase3xNs = positive_mod(relativeNs - phaseOriginNs, period3xNs);
        fill_profile(folded3xProfiles, folded3xTotals, channel, phase3xNs, folded3xBins, 0.0, period3xNs);
        if (channel == config.scoreChannel) {
            fill_phase_tot(folded3xPhaseTotHistograms, channel, phase3xNs, folded3xBins, 0.0, period3xNs, totNs, config.totBins, config.totMinNs, config.totMaxNs);
            fill_sparse_histogram(
                folded3xPhaseCh0TimeHistograms,
                channel,
                phase3xNs,
                folded3xBins,
                0.0,
                period3xNs,
                relativeNs,
                config.ch0RefBins,
                config.ch0RefMinNs,
                config.ch0RefMaxNs
            );
            fill_sparse_histogram(
                folded3xPhaseTriggerTimeHistograms,
                channel,
                phase3xNs,
                folded3xBins,
                0.0,
                period3xNs,
                triggerRelativeNs,
                config.triggerRefBins,
                config.triggerRefMinNs,
                config.triggerRefMaxNs
            );
        }
        stats.folded3xAcceptedPulses += 1;
        stats.folded3xAcceptedByChannel[static_cast<std::size_t>(channel)] += 1;
    };

    fill_all_modes(0, 0.0, referenceHit->riseNs, referenceHit->totNs);
    for (int channel = 1; channel < kChannelCount; ++channel) {
        for (const PulseHit& hit : window.hitsByChannel[static_cast<std::size_t>(channel)]) {
            fill_all_modes(channel, hit.riseNs - referenceHit->riseNs, hit.riseNs, hit.totNs);
        }
    }
}

void write_profile_matrix(const fs::path& outputPath,
                         const std::string& referenceLabel,
                         double min,
                         double max,
                         int bins,
                         const Config& config,
                         const Stats& stats,
                         const CandidateMetrics& best,
                         const ChannelProfiles& profiles,
                         const std::vector<std::uint64_t>& totals,
                         int foldFactor) {
    std::ofstream output(outputPath);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write RF profile matrix: " + outputPath.string());
    }

    output << "# input_pulse_table=" << config.inputPath << '\n';
    if (!stats.sourceInputFile.empty()) {
        output << "# input_raw_file=" << stats.sourceInputFile << '\n';
    }
    output << "# run_key=" << config.runKey << '\n';
    output << "# x_reference=" << referenceLabel << '\n';
    output << "# x_min_ns=" << format_double(min) << '\n';
    output << "# x_max_ns=" << format_double(max) << '\n';
    output << "# bins=" << bins << '\n';
    output << "# deduced_period_ns=" << format_double(best.periodNs) << '\n';
    output << "# phase_origin_ns=" << format_double(best.phaseOriginNs) << '\n';
    output << "# peak_center_ns=" << format_double(best.peakCenterNs) << '\n';
    output << "# fold_factor=" << foldFactor << '\n';
    output << "# score_channel=" << config.scoreChannel << '\n';
    output << "# score_tot_min_ns=" << format_double(config.scoreTotMinNs) << '\n';
    output << "# score_tot_max_ns=" << format_double(config.scoreTotMaxNs) << '\n';
    output << "# ch0_valid_rise_min_ns=" << format_double(config.ch0ValidRiseMinNs) << '\n';
    output << "# ch0_valid_rise_max_ns=" << format_double(config.ch0ValidRiseMaxNs) << '\n';
    output << "# ch0_valid_tot_min_ns=" << format_double(config.ch0ValidTotMinNs) << '\n';
    output << "# ch0_valid_tot_max_ns=" << format_double(config.ch0ValidTotMaxNs) << '\n';
    output << "# valid_ch0_windows=" << stats.validCh0Windows << '\n';
    output << "# windows_without_valid_ch0_candidate=" << stats.windowsWithoutValidCh0Candidate << '\n';
    output << "# columns: bin_index time_center_ns";
    for (int channel = 0; channel < kChannelCount; ++channel) {
        output << " ch" << std::setfill('0') << std::setw(2) << channel;
    }
    output << " total\n";
    output << std::setfill(' ');

    for (int index = 0; index < bins; ++index) {
        const double center = min + (static_cast<double>(index) + 0.5) * ((max - min) / static_cast<double>(bins));
        output << index << '\t' << format_double(center);
        for (int channel = 0; channel < kChannelCount; ++channel) {
            output << '\t' << profiles[static_cast<std::size_t>(channel)][static_cast<std::size_t>(index)];
        }
        output << '\t' << totals[static_cast<std::size_t>(index)] << '\n';
    }
}

void write_phase_tot_sparse(const fs::path& outputPath,
                            const std::string& xReference,
                            double xMin,
                            double xMax,
                            int xBins,
                            const Config& config,
                            const Stats& stats,
                            const CandidateMetrics& best,
                            const ChannelPhaseTotHistograms& histograms,
                            int foldFactor) {
    std::ofstream output(outputPath);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write folded phase-ToT histogram: " + outputPath.string());
    }

    output << "# input_pulse_table=" << config.inputPath << '\n';
    if (!stats.sourceInputFile.empty()) {
        output << "# input_raw_file=" << stats.sourceInputFile << '\n';
    }
    output << "# run_key=" << config.runKey << '\n';
    output << "# x_reference=" << xReference << '\n';
    output << "# x_min_ns=" << format_double(xMin) << '\n';
    output << "# x_max_ns=" << format_double(xMax) << '\n';
    output << "# x_bins=" << xBins << '\n';
    output << "# y_reference=ToT (ns)\n";
    output << "# y_min_ns=" << format_double(config.totMinNs) << '\n';
    output << "# y_max_ns=" << format_double(config.totMaxNs) << '\n';
    output << "# y_bins=" << config.totBins << '\n';
    output << "# tot_min_ns=" << format_double(config.totMinNs) << '\n';
    output << "# tot_max_ns=" << format_double(config.totMaxNs) << '\n';
    output << "# tot_bins=" << config.totBins << '\n';
    output << "# channel_count=" << kChannelCount << '\n';
    output << "# deduced_period_ns=" << format_double(best.periodNs) << '\n';
    output << "# phase_origin_ns=" << format_double(best.phaseOriginNs) << '\n';
    output << "# peak_center_ns=" << format_double(best.peakCenterNs) << '\n';
    output << "# fold_factor=" << foldFactor << '\n';
    output << "# columns: channel x_bin y_bin count\n";

    for (int channel = 0; channel < kChannelCount; ++channel) {
        const auto& storage = histograms[static_cast<std::size_t>(channel)];
        if (storage.empty()) {
            continue;
        }
        for (int yBin = 0; yBin < config.totBins; ++yBin) {
            for (int xBin = 0; xBin < xBins; ++xBin) {
                const std::size_t flatIndex = static_cast<std::size_t>(yBin) * static_cast<std::size_t>(xBins)
                                            + static_cast<std::size_t>(xBin);
                const std::uint64_t count = storage[flatIndex];
                if (count == 0) {
                    continue;
                }
                output << channel << ' ' << xBin << ' ' << yBin << ' ' << count << '\n';
            }
        }
    }
}

void write_phase_value_sparse(const fs::path& outputPath,
                              const std::string& xReference,
                              const std::string& yReference,
                              double xMin,
                              double xMax,
                              int xBins,
                              double yMin,
                              double yMax,
                              int yBins,
                              const Config& config,
                              const Stats& stats,
                              const CandidateMetrics& best,
                              const SparseHistogram& histogram,
                              int foldFactor) {
    std::ofstream output(outputPath);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write folded sparse histogram: " + outputPath.string());
    }

    output << "# input_pulse_table=" << config.inputPath << '\n';
    if (!stats.sourceInputFile.empty()) {
        output << "# input_raw_file=" << stats.sourceInputFile << '\n';
    }
    output << "# run_key=" << config.runKey << '\n';
    output << "# x_reference=" << xReference << '\n';
    output << "# y_reference=" << yReference << '\n';
    output << "# x_min_ns=" << format_double(xMin) << '\n';
    output << "# x_max_ns=" << format_double(xMax) << '\n';
    output << "# x_bins=" << xBins << '\n';
    output << "# y_min_ns=" << format_double(yMin) << '\n';
    output << "# y_max_ns=" << format_double(yMax) << '\n';
    output << "# y_bins=" << yBins << '\n';
    output << "# channel_count=" << kChannelCount << '\n';
    output << "# deduced_period_ns=" << format_double(best.periodNs) << '\n';
    output << "# phase_origin_ns=" << format_double(best.phaseOriginNs) << '\n';
    output << "# peak_center_ns=" << format_double(best.peakCenterNs) << '\n';
    output << "# fold_factor=" << foldFactor << '\n';
    output << "# columns: channel x_bin_index_zero_based y_bin_index_zero_based count\n";

    for (int channel = 0; channel < kChannelCount; ++channel) {
        const auto& storage = histogram[static_cast<std::size_t>(channel)];
        if (storage.empty()) {
            continue;
        }
        for (int xBin = 0; xBin < xBins; ++xBin) {
            for (int yBin = 0; yBin < yBins; ++yBin) {
                const std::size_t key = static_cast<std::size_t>(sparse_key(xBin, yBin, yBins));
                const std::uint32_t count = storage[key];
                if (count == 0U) {
                    continue;
                }
                output << channel << ' ' << xBin << ' ' << yBin << ' ' << count << '\n';
            }
        }
    }
}

void write_phase_segments_file(const fs::path& outputPath,
                               const Config& config,
                               const CandidateMetrics& best,
                               const std::vector<PhaseSegment>& segments) {
    std::ofstream output(outputPath);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write RF phase segments file: " + outputPath.string());
    }

    output << "# input_pulse_table=" << config.inputPath << '\n';
    output << "# run_key=" << config.runKey << '\n';
    output << "# deduced_period_ns=" << format_double(best.periodNs) << '\n';
    output << "# global_phase_origin_ns=" << format_double(best.phaseOriginNs) << '\n';
    output << "# segment_block_score_pulses=" << config.segmentBlockScorePulses << '\n';
    output << "# segment_persistent_blocks=" << config.segmentPersistentBlocks << '\n';
    output << "# segment_phase_jump_ns=" << format_double(config.segmentPhaseJumpNs) << '\n';
    output << "# segment_gap_factor=" << format_double(config.segmentGapFactor) << '\n';
    output << "# segment_min_gap_seconds=" << format_double(config.segmentMinGapSeconds) << '\n';
    output << "# columns: segment_id start_global_time_seconds end_global_time_seconds phase_origin_ns score_pulses phase_shift_from_previous_ns reason\n";
    for (const PhaseSegment& segment : segments) {
        output << segment.segmentId << ' '
               << format_double(segment.startGlobalTimeSeconds) << ' '
               << format_double(segment.endGlobalTimeSeconds) << ' '
               << format_double(segment.phaseOriginNs) << ' '
               << segment.scorePulses << ' '
               << format_double(segment.phaseShiftFromPreviousNs) << ' '
               << segment.reason << '\n';
    }
}

void write_summary(const fs::path& outputPath,
                   const Config& config,
                   const Stats& stats,
                   const CandidateMetrics& best,
                   const std::vector<PhaseSegment>& phaseSegments,
                   const fs::path& ch0RefPath,
                   const fs::path& foldedPath,
                   const fs::path& folded3xPath,
                   const fs::path& foldedPhaseTotPath,
                   const fs::path& folded3xPhaseTotPath,
                   const fs::path& foldedPhaseCh0TimePath,
                   const fs::path& foldedPhaseTriggerTimePath,
                   const fs::path& folded3xPhaseCh0TimePath,
                   const fs::path& folded3xPhaseTriggerTimePath,
                   const fs::path& phaseSegmentsPath,
                   const fs::path& scanPrefix) {
    std::ofstream output(outputPath);
    if (!output.is_open()) {
        throw std::runtime_error("Unable to write cleaned RF summary: " + outputPath.string());
    }

    output << "input_pulse_table=" << config.inputPath << '\n';
    if (!stats.sourceInputFile.empty()) {
        output << "input_raw_file=" << stats.sourceInputFile << '\n';
    }
    output << "run_key=" << config.runKey << '\n';
    output << "ch0_ref_output=" << ch0RefPath << '\n';
    output << "folded_output=" << foldedPath << '\n';
    output << "folded3x_output=" << folded3xPath << '\n';
    output << "folded_phase_tot_output=" << foldedPhaseTotPath << '\n';
    output << "folded3x_phase_tot_output=" << folded3xPhaseTotPath << '\n';
    output << "folded_phase_ch0_ref_time_output=" << foldedPhaseCh0TimePath << '\n';
    output << "folded_phase_trigger_time_output=" << foldedPhaseTriggerTimePath << '\n';
    output << "folded3x_phase_ch0_ref_time_output=" << folded3xPhaseCh0TimePath << '\n';
    output << "folded3x_phase_trigger_time_output=" << folded3xPhaseTriggerTimePath << '\n';
    output << "rf_segment_output=" << phaseSegmentsPath << '\n';
    output << "scan_prefix=" << scanPrefix << '\n';
    output << "score_channel=" << config.scoreChannel << '\n';
    output << "score_tot_min_ns=" << format_double(config.scoreTotMinNs) << '\n';
    output << "score_tot_max_ns=" << format_double(config.scoreTotMaxNs) << '\n';
    output << "deduced_period_ns=" << format_double(best.periodNs) << '\n';
    output << "phase_origin_ns=" << format_double(best.phaseOriginNs) << '\n';
    output << "rf_phase_mode=" << (phaseSegments.size() > 1 ? "segmented" : "global") << '\n';
    output << "rf_segment_count=" << phaseSegments.size() << '\n';
    output << "rf_segment_guard_triggered=" << (phaseSegments.size() > 1 ? "true" : "false") << '\n';
    output << "rf_segment_block_score_pulses=" << config.segmentBlockScorePulses << '\n';
    output << "rf_segment_persistent_blocks=" << config.segmentPersistentBlocks << '\n';
    output << "rf_segment_phase_jump_ns=" << format_double(config.segmentPhaseJumpNs) << '\n';
    output << "rf_segment_gap_factor=" << format_double(config.segmentGapFactor) << '\n';
    output << "rf_segment_min_gap_seconds=" << format_double(config.segmentMinGapSeconds) << '\n';
    output << "single_rf_segment=" << (config.singleRfSegment ? "true" : "false") << '\n';
    output << "skip_cycle_residual_diagnostics=" << (config.skipCycleResidualDiagnostics ? "true" : "false") << '\n';
    output << "peak_center_ns=" << format_double(best.peakCenterNs) << '\n';
    output << "peak_height=" << best.peakHeight << '\n';
    output << "selected_pulses=" << best.selectedPulses << '\n';
    output << "selected_fraction=" << format_double(best.selectedFraction) << '\n';
    output << "sigma_ns=" << format_double(best.sigmaNs) << '\n';
    output << "coherence=" << format_double(best.coherence) << '\n';
    output << "drift_slope_ns_per_cycle=" << format_double(best.driftSlopeNsPerCycle) << '\n';
    output << "merit=" << format_double(best.merit) << '\n';
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
    output << "score_pulses_before_stride=" << stats.scorePulsesBeforeStride << '\n';
    output << "score_pulses_stored=" << stats.scorePulsesStored << '\n';
    output << "ch0_ref_accepted_pulses=" << stats.ch0RefAcceptedPulses << '\n';
    output << "folded_accepted_pulses=" << stats.foldedAcceptedPulses << '\n';
    output << "folded3x_accepted_pulses=" << stats.folded3xAcceptedPulses << '\n';
    output << "ch0_ref_rejected_outside_range=" << stats.ch0RefRejectedOutsideRange << '\n';
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
        std::vector<ScoredPulse> scorePulses;
        scorePulses.reserve(1 << 20);
        WindowState currentWindow;
        std::uint64_t scoreOrdinal = 0;

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

            if (!currentWindow.hasValue) {
                currentWindow.hasValue = true;
                currentWindow.windowIndex = row.windowIndex;
                currentWindow.hasGlobalTriggerSeconds = true;
                currentWindow.globalTriggerSeconds = row.globalTriggerSeconds;
            }
            if (row.windowIndex != currentWindow.windowIndex) {
                finalize_window_for_scan(currentWindow, config, stats, scorePulses, scoreOrdinal);
                currentWindow = WindowState{};
                currentWindow.hasValue = true;
                currentWindow.windowIndex = row.windowIndex;
                currentWindow.hasGlobalTriggerSeconds = true;
                currentWindow.globalTriggerSeconds = row.globalTriggerSeconds;
            }
            currentWindow.hitsByChannel[static_cast<std::size_t>(row.channel)].push_back(PulseHit{row.riseNs, row.totNs});
        }
        finalize_window_for_scan(currentWindow, config, stats, scorePulses, scoreOrdinal);

        if (config.runKey.empty()) {
            config.runKey = strip_known_suffixes(fs::path(config.inputPath).filename().string());
        }

        const std::vector<CandidateMetrics> scanRows = scan_candidates(config, scorePulses);
        const CandidateMetrics best = best_valid_candidate(scanRows);
        const Histogram1D bestPhaseProfile = build_best_phase_profile(config, scorePulses, best);
        Histogram2D bestPhaseTot = make_best_phase_tot_histogram(config, best);
        const std::vector<CycleResidualStats> cycleResiduals = config.skipCycleResidualDiagnostics
            ? std::vector<CycleResidualStats>{}
            : build_cycle_residuals(config, scorePulses, best);
        const std::vector<CycleResidualPoint> cycleResidualPoints = config.skipCycleResidualDiagnostics
            ? std::vector<CycleResidualPoint>{}
            : build_cycle_residual_points(config, scorePulses, best);
        const std::vector<PhaseBlock> phaseBlocks = config.singleRfSegment
            ? std::vector<PhaseBlock>{}
            : build_phase_blocks(config, scorePulses, best.periodNs);
        const std::vector<PhaseSegment> phaseSegments = config.singleRfSegment
            ? std::vector<PhaseSegment>{build_single_phase_segment(scorePulses, best)}
            : build_phase_segments(config, scorePulses, phaseBlocks, best);

        const int foldedBins = std::max(16, static_cast<int>(std::llround(best.periodNs / config.phaseBinWidthNs)));
        const int folded3xBins = 3 * foldedBins;
        ChannelProfiles ch0RefProfiles;
        ChannelProfiles foldedProfiles;
        ChannelProfiles folded3xProfiles;
        ChannelPhaseTotHistograms foldedPhaseTotHistograms;
        SparseHistogram foldedPhaseCh0TimeHistograms;
        SparseHistogram foldedPhaseTriggerTimeHistograms;
        ChannelPhaseTotHistograms folded3xPhaseTotHistograms;
        SparseHistogram folded3xPhaseCh0TimeHistograms;
        SparseHistogram folded3xPhaseTriggerTimeHistograms;
        ensure_profile_sizes(ch0RefProfiles, config.ch0RefBins);
        ensure_profile_sizes(foldedProfiles, foldedBins);
        ensure_profile_sizes(folded3xProfiles, folded3xBins);
        ensure_phase_tot_sizes(foldedPhaseTotHistograms, foldedBins, config.totBins, config.scoreChannel);
        ensure_phase_tot_sizes(folded3xPhaseTotHistograms, folded3xBins, config.totBins, config.scoreChannel);
        ensure_sparse_histogram_sizes(foldedPhaseCh0TimeHistograms, foldedBins, config.ch0RefBins, config.scoreChannel);
        ensure_sparse_histogram_sizes(foldedPhaseTriggerTimeHistograms, foldedBins, config.triggerRefBins, config.scoreChannel);
        ensure_sparse_histogram_sizes(folded3xPhaseCh0TimeHistograms, folded3xBins, config.ch0RefBins, config.scoreChannel);
        ensure_sparse_histogram_sizes(folded3xPhaseTriggerTimeHistograms, folded3xBins, config.triggerRefBins, config.scoreChannel);
        std::vector<std::uint64_t> ch0RefTotals(static_cast<std::size_t>(config.ch0RefBins), 0);
        std::vector<std::uint64_t> foldedTotals(static_cast<std::size_t>(foldedBins), 0);
        std::vector<std::uint64_t> folded3xTotals(static_cast<std::size_t>(folded3xBins), 0);

        input.close();
        drop_file_cache_hint(config.inputPath);

        std::ifstream profileInput(config.inputPath);
        if (!profileInput.is_open()) {
            throw std::runtime_error("Unable to reopen cleaned pulse table for profile pass: " + config.inputPath);
        }

        currentWindow = WindowState{};
        while (std::getline(profileInput, line)) {
            if (line.empty()) {
                continue;
            }
            std::string key;
            std::string value;
            if (parse_metadata_line(line, key, value) || starts_with(line, "window_index")) {
                continue;
            }
            PulseRow row;
            if (!parse_pulse_row(line, row)) {
                continue;
            }

            if (!currentWindow.hasValue) {
                currentWindow.hasValue = true;
                currentWindow.windowIndex = row.windowIndex;
                currentWindow.hasGlobalTriggerSeconds = true;
                currentWindow.globalTriggerSeconds = row.globalTriggerSeconds;
            }
            if (row.windowIndex != currentWindow.windowIndex) {
                finalize_window_for_profiles(currentWindow,
                                             config,
                                             best,
                                             bestPhaseTot,
                                             phaseSegments,
                                             stats,
                                             ch0RefProfiles,
                                             ch0RefTotals,
                                             foldedProfiles,
                                             foldedTotals,
                                             foldedPhaseTotHistograms,
                                             foldedPhaseCh0TimeHistograms,
                                             foldedPhaseTriggerTimeHistograms,
                                             folded3xProfiles,
                                             folded3xTotals,
                                             folded3xPhaseTotHistograms,
                                             folded3xPhaseCh0TimeHistograms,
                                             folded3xPhaseTriggerTimeHistograms);
                currentWindow = WindowState{};
                currentWindow.hasValue = true;
                currentWindow.windowIndex = row.windowIndex;
                currentWindow.hasGlobalTriggerSeconds = true;
                currentWindow.globalTriggerSeconds = row.globalTriggerSeconds;
            }
            currentWindow.hitsByChannel[static_cast<std::size_t>(row.channel)].push_back(PulseHit{row.riseNs, row.totNs});
        }
        finalize_window_for_profiles(currentWindow,
                                     config,
                                     best,
                                     bestPhaseTot,
                                     phaseSegments,
                                     stats,
                                     ch0RefProfiles,
                                     ch0RefTotals,
                                     foldedProfiles,
                                     foldedTotals,
                                     foldedPhaseTotHistograms,
                                     foldedPhaseCh0TimeHistograms,
                                     foldedPhaseTriggerTimeHistograms,
                                     folded3xProfiles,
                                     folded3xTotals,
                                     folded3xPhaseTotHistograms,
                                     folded3xPhaseCh0TimeHistograms,
                                     folded3xPhaseTriggerTimeHistograms);
                        profileInput.close();
                        drop_file_cache_hint(config.inputPath);

        const fs::path runDir = fs::path(config.outputRoot) / config.runKey;
        const fs::path ch0RefDir = runDir / config.ch0RefDirName;
        const fs::path foldedDir = runDir / config.foldedDirName;
        const fs::path folded3xDir = runDir / config.folded3xDirName;
        const fs::path scanDir = runDir / config.scanDirName;
        fs::create_directories(ch0RefDir);
        fs::create_directories(foldedDir);
        fs::create_directories(folded3xDir);
        fs::create_directories(scanDir);

        const fs::path ch0RefPath = ch0RefDir / (config.runKey + "_Ch0_ref_Rates_scan_matrix.tsv");
        const fs::path foldedPath = foldedDir / (config.runKey + "_Folded_RF_matrix.tsv");
        const fs::path folded3xPath = folded3xDir / (config.runKey + "_Folded_RF_3x_matrix.tsv");
        const fs::path foldedPhaseTotPath = foldedDir / (config.runKey + "_Folded_RF_phase_tot_sparse.tsv");
        const fs::path foldedPhaseCh0TimePath = foldedDir / (config.runKey + "_Folded_RF_phase_vs_ch0_ref_time_sparse.tsv");
        const fs::path foldedPhaseTriggerTimePath = foldedDir / (config.runKey + "_Folded_RF_phase_vs_trigger_time_sparse.tsv");
        const fs::path folded3xPhaseTotPath = folded3xDir / (config.runKey + "_Folded_RF_3x_phase_tot_sparse.tsv");
        const fs::path folded3xPhaseCh0TimePath = folded3xDir / (config.runKey + "_Folded_RF_3x_phase_vs_ch0_ref_time_sparse.tsv");
        const fs::path folded3xPhaseTriggerTimePath = folded3xDir / (config.runKey + "_Folded_RF_3x_phase_vs_trigger_time_sparse.tsv");
        const fs::path scanPrefix = scanDir / (config.runKey + "_rf_period_scan");
        const fs::path scanPath = scanPrefix.string() + "_scan.txt";
        const fs::path bestPhaseProfilePath = scanPrefix.string() + "_best_phase_profile_hist.txt";
        const fs::path bestPhaseTotPath = scanPrefix.string() + "_best_phase_tot_hist.txt";
        const fs::path cycleResidualPath = scanPrefix.string() + "_best_cycle_residuals.txt";
        const fs::path cycleResidualPointsPath = scanPrefix.string() + "_best_cycle_residual_points.txt";
        const fs::path phaseSegmentsPath = scanPrefix.string() + "_phase_segments.tsv";
        const fs::path summaryPath = runDir / (config.runKey + "_cleaned_rf_period_summary.txt");

        write_profile_matrix(ch0RefPath,
                             "validated_first_ch0_reference_hit",
                             config.ch0RefMinNs,
                             config.ch0RefMaxNs,
                             config.ch0RefBins,
                             config,
                             stats,
                             best,
                             ch0RefProfiles,
                             ch0RefTotals,
                             0);
        write_profile_matrix(foldedPath,
                             "deduced_rf_period_folded",
                             0.0,
                             best.periodNs,
                             foldedBins,
                             config,
                             stats,
                             best,
                             foldedProfiles,
                             foldedTotals,
                             1);
        write_profile_matrix(folded3xPath,
                             "deduced_rf_period_folded_3x",
                             0.0,
                             3.0 * best.periodNs,
                             folded3xBins,
                             config,
                             stats,
                             best,
                             folded3xProfiles,
                             folded3xTotals,
                             3);
                    write_phase_tot_sparse(foldedPhaseTotPath,
                                   "deduced_rf_period_folded",
                                   0.0,
                                   best.periodNs,
                                   foldedBins,
                                   config,
                                   stats,
                                   best,
                                   foldedPhaseTotHistograms,
                                   1);
                    write_phase_value_sparse(foldedPhaseCh0TimePath,
                                   "deduced_rf_period_folded",
                                   "time relative to validated first ch0 reference hit",
                                   0.0,
                                   best.periodNs,
                                   foldedBins,
                                   config.ch0RefMinNs,
                                   config.ch0RefMaxNs,
                                   config.ch0RefBins,
                                   config,
                                   stats,
                                   best,
                                   foldedPhaseCh0TimeHistograms,
                                   1);
                    write_phase_value_sparse(foldedPhaseTriggerTimePath,
                                   "deduced_rf_period_folded",
                                   "time relative to trigger window",
                                   0.0,
                                   best.periodNs,
                                   foldedBins,
                                   config.triggerRefMinNs,
                                   config.triggerRefMaxNs,
                                   config.triggerRefBins,
                                   config,
                                   stats,
                                   best,
                                   foldedPhaseTriggerTimeHistograms,
                                   1);
                    write_phase_tot_sparse(folded3xPhaseTotPath,
                                   "deduced_rf_period_folded_3x",
                                   0.0,
                                   3.0 * best.periodNs,
                                   folded3xBins,
                                   config,
                                   stats,
                                   best,
                                   folded3xPhaseTotHistograms,
                                   3);
                    write_phase_value_sparse(folded3xPhaseCh0TimePath,
                                   "deduced_rf_period_folded_3x",
                                   "time relative to validated first ch0 reference hit",
                                   0.0,
                                   3.0 * best.periodNs,
                                   folded3xBins,
                                   config.ch0RefMinNs,
                                   config.ch0RefMaxNs,
                                   config.ch0RefBins,
                                   config,
                                   stats,
                                   best,
                                   folded3xPhaseCh0TimeHistograms,
                                   3);
                    write_phase_value_sparse(folded3xPhaseTriggerTimePath,
                                   "deduced_rf_period_folded_3x",
                                   "time relative to trigger window",
                                   0.0,
                                   3.0 * best.periodNs,
                                   folded3xBins,
                                   config.triggerRefMinNs,
                                   config.triggerRefMaxNs,
                                   config.triggerRefBins,
                                   config,
                                   stats,
                                   best,
                                   folded3xPhaseTriggerTimeHistograms,
                                   3);
                    write_scan_file(scanPath, config, scanRows);
                    write_histogram_1d(bestPhaseProfilePath, bestPhaseProfile);
                    write_histogram_2d(bestPhaseTotPath, bestPhaseTot);
                    write_cycle_residual_file(cycleResidualPath, best, cycleResiduals);
                    write_cycle_residual_points_file(cycleResidualPointsPath, best, cycleResidualPoints);
                    write_phase_segments_file(phaseSegmentsPath, config, best, phaseSegments);
                    write_summary(summaryPath,
                                  config,
                                  stats,
                                  best,
                              phaseSegments,
                                  ch0RefPath,
                                  foldedPath,
                                  folded3xPath,
                                  foldedPhaseTotPath,
                                  folded3xPhaseTotPath,
                                  foldedPhaseCh0TimePath,
                                  foldedPhaseTriggerTimePath,
                                  folded3xPhaseCh0TimePath,
                                  folded3xPhaseTriggerTimePath,
                              phaseSegmentsPath,
                                  scanPrefix);

        std::cout << "Ch0-ref matrix: " << ch0RefPath << '\n';
        std::cout << "Folded matrix: " << foldedPath << '\n';
        std::cout << "Folded 3x matrix: " << folded3xPath << '\n';
                    std::cout << "Folded phase-ToT sparse: " << foldedPhaseTotPath << '\n';
                    std::cout << "Folded phase-vs-ch0-ref-time sparse: " << foldedPhaseCh0TimePath << '\n';
                    std::cout << "Folded phase-vs-trigger-time sparse: " << foldedPhaseTriggerTimePath << '\n';
                    std::cout << "Folded 3x phase-ToT sparse: " << folded3xPhaseTotPath << '\n';
                    std::cout << "Folded 3x phase-vs-ch0-ref-time sparse: " << folded3xPhaseCh0TimePath << '\n';
                    std::cout << "Folded 3x phase-vs-trigger-time sparse: " << folded3xPhaseTriggerTimePath << '\n';
                    std::cout << "Scan prefix: " << scanPrefix << '\n';
        std::cout << "Summary: " << summaryPath << '\n';
        std::cout << "Best period (ns): " << format_double(best.periodNs) << '\n';
        std::cout << "Valid ch0 windows: " << stats.validCh0Windows << '\n';
        std::cout << "Windows without valid ch0 candidate: " << stats.windowsWithoutValidCh0Candidate << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
