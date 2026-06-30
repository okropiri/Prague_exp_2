# NCAL RF Phase Scan Walkthrough

## Why This Exists

This note explains the RF-period search in plain language.

It is meant to answer one question:

How do we go from raw `ch0`-referenced `NCAL` times to one best trial period such as about `38.85 ns`?

## What Happens Before The Search

1. Read the raw DOGMA file.
2. Keep only the signal windows.
3. In each signal window, take the earliest `ch0` rise as the local time zero.
4. Pair each `NCAL` rise with its matching fall to get a `ToT`.
5. Compute the raw `ch0`-referenced `NCAL` time:

   $$
   t_{\mathrm{rel}} = t_{\mathrm{NCAL\ rise}} - t_{\mathrm{earliest\ ch0}}
   $$

6. Keep only pulses that pass the accepted time and `ToT` cuts.
7. Store all accepted pulses for the final outputs.
8. For the search only, keep every `N`th accepted pulse if `scoreStride > 1`.

Important:

- The search may use a subsample for speed.
- The final phase-vs-ToT plot is built later from all accepted stored pulses.

## Current Search Knobs

These are the main search settings used by the current tool.

- `initialPeriodMinNs = 39.0`
- `initialPeriodMaxNs = 41.0`
- `initialStepNs = 0.05`
- `refineRounds = 3`
- `refineHalfSpanSteps = 5`
- `refineFactor = 10.0`
- `phaseBinWidthNs = 0.25`
- `peakWindowNs = 8.0`
- `minSelectedPulses = 1000`
- recentering passes per trial period = `2`
- `scoreStride = 1` by default

For the first `3.4 m / 0000` run, the search used:

- `scoreStride = 32`
- `peakWindowNs = 8.0`, which means a full selection width of `16 ns`

## What The Code Does For One Trial Period

Take one tested period, call it `P`.

### Step 1: Fold The Sampled Times

For every sampled raw `NCAL` time `t`, compute a folded phase:

$$
\phi = t \bmod P
$$

Now every sampled pulse lies between `0` and `P`.

### Step 2: Build A 1D Phase Histogram

Ignore `ToT` for the moment.

Take all folded phases and count how many fall into each small phase bin.

This gives a temporary 1D histogram: phase counts versus phase.

### Step 3: Search The Full Phase Range For Seed Windows

Do not lock onto only one tallest phase bin.

Instead, scan seed centers across the whole folded phase.

For each seed center, count how many sampled pulses lie inside the full `16 ns` window:

$$
|r| \le 8 \text{ ns}
$$

This is a circular search over the whole `0` to `P` folded phase, so a broader real blob is not discarded just because some narrower feature owns the single tallest bin.

The code keeps a small set of the strongest non-overlapping seed windows and tests each one.

So at this point, the search has already looked across the full phase range before any local recentering starts.

### Step 4: Recenter The Peak Twice

Start from one tested seed-window center.

Then do this two times:

1. Go through all sampled folded phases again.
2. For each pulse, compute its signed wrapped residual from the current center.
3. Keep only pulses inside the peak window:

   $$
   |r| \le 8 \text{ ns}
   $$

4. Average those signed wrapped residuals.
5. Shift the current peak center by that average.

This recentering step is trying to make the mean signed wrapped residual close to zero.

Important:

- the recentering itself is still local,
- but the initial seed search is no longer tied to only the tallest single histogram bin.

What that residual means:

- First, fold the event time into one trial period.
- Then compare it with the current peak center.
- The result is a signed distance from the center of the peak.

So the residual is not the raw `NCAL - ch0` time itself.
It is the small offset of that `ch0`-referenced event from the expected phase center for the tested period.

Important:

- This uses signed wrapped residuals, not absolute distances.
- It is a short recentering step, not a long convergence loop.

### Short Picture Of A Centered Residual

Suppose the trial period is about `38.85 ns`.

Then the expected repeated phase centers are separated by about:

$$
..., -38.85,
0,
38.85,
77.70,
...
$$

Now take one `ch0`-referenced event time.

- If it lands exactly on its expected phase center, the residual is near `0 ns`.
- If it is a bit early, the residual is negative.
- If it is a bit late, the residual is positive.

That is why the code converts the phase difference into a centered wrapped residual:

$$
r \in [-P/2, +P/2]
$$

This centered view is much easier to interpret than a raw phase between `0` and `P`.

Example:

- event near `38.10 ns` for a center at `38.85 ns` gives about `-0.75 ns`
- event near `39.20 ns` for a center at `38.85 ns` gives about `+0.35 ns`

So the residual always means:

How far is this event from the nearest expected RF-cycle center?

### Step 5: Make The Final Selection For This Trial Period

After the center has been recentered, go through the sampled pulses again.

Keep only pulses whose wrapped residual still satisfies:

$$
|r| \le 8 \text{ ns}
$$

This is the selected cluster for this trial period.

At this stage:

- `Nscore` means all sampled pulses that entered the test of this trial period
- `Nselected` means the sampled pulses that fall inside the final `\pm 8 ns` peak window

### Step 6: Compute The Trial-Period Quality Numbers

From that selected cluster, compute:

- `selectedPulses = Nselected`: how many sampled pulses are inside the peak window
- `selectedFraction`:

  $$
  \text{selectedFraction} = \frac{N_{\mathrm{selected}}}{N_{\mathrm{score}}}
  $$

- `meanResidualNs`: the mean signed wrapped residual of the selected pulses

   $$
   \mu = \frac{1}{N_{\mathrm{selected}}} \sum r_i
   $$

- `sigmaNs`: the standard deviation of the selected wrapped residuals

   $$
   \sigma^2 = \frac{1}{N_{\mathrm{selected}}} \sum r_i^2 - \mu^2
   $$

   $$
   \sigma = \sqrt{\sigma^2}
   $$

- `driftSlopeNsPerCycle`: whether residuals walk with cycle index

To define that slope, the code also assigns each selected event a cycle index.

The cycle index is the relative RF-cycle label of that event for the current tested period:

- cycle `0` means the cycle closest to the chosen phase center
- cycle `1` means one period later
- cycle `-1` means one period earlier

So the cycle index is not an absolute cycle count over the whole run.
It is a relative cycle number inside the `ch0`-referenced analysis window.

For the current settings, that window is:

$$
[-6000 \text{ ns}, +6000 \text{ ns})
$$

which is a total width of `12000 ns`.

With a period near `38.85 ns`, this means only a few hundred RF cycles fit into the kept time range:

$$
12000 / 38.85 \approx 309
$$

So when a plot shows about a few hundred cycle indices, that is expected.
It is showing the cycle slots inside the selected trigger window, not the total number of cycles in the whole acquisition.

Then compute the main ranking score:

$$
\text{merit} = \frac{\text{selectedFraction}}{\sigma}
$$

Because `Nscore` is fixed during one scan, this is effectively proportional to:

$$
\frac{N_{\mathrm{selected}}}{\sigma}
$$

So the preferred trial period is one that keeps many pulses inside one peak and makes that peak narrow.

### What The Drift Slope Means

The drift slope is the fitted change of residual per RF cycle:

$$
a = \frac{\Delta r}{\Delta n}
$$

Its unit is:

$$
\mathrm{ns/cycle}
$$

Interpretation:

- slope near `0`: the selected peak stays centered as cycle index changes
- positive slope: later cycles tend to shift to more positive residual
- negative slope: later cycles tend to shift to more negative residual

So this quantity measures whether the selected cluster drifts across the analysis window.

If the tested period is slightly wrong, the phase error accumulates from cycle to cycle.
Then the residual tends to change linearly with cycle index, and the fitted slope moves away from zero.

If the tested period is close to the true one, that accumulated drift becomes small, and the slope tends to move back toward zero.

This is why the slope panel in the period scan often changes sign near the best period.

### Step 7: Store One Scan Point

That one trial period contributes one point to the scan outputs:

- period
- merit
- sigma
- drift slope
- selected fraction

Then the tool moves to the next tested period.

## What The Multi-Round Scan Means

The scan is not done only once.

It runs in rounds:

1. coarse scan over a broad period range
2. first refinement around the best coarse point
3. second refinement around the best refined point

So the later rounds zoom in on the most promising region.

## How To Read The Main Period-Scan Panels

The scan plot has four stacked panels.

### Selected Fraction

This shows how many sampled pulses survive the final peak-window cut for each tested period.

Higher is better because more events join one consistent folded peak.

### Sigma

This shows the width of the selected residual cluster.

Lower is better because the folded peak is tighter.

### Merit

This is the ranking score:

$$
	ext{merit} = \frac{\text{selectedFraction}}{\sigma}
$$

So the merit prefers trial periods that keep many events while also making the selected cluster narrow.

### Drift Slope

This is the fitted residual-versus-cycle slope for the same selected cluster.

Important:

- it is not the derivative of merit
- it is not computed from the merit curve

The merit and slope can still look related because both respond to the same underlying effect:

- if the trial period is wrong, phase error accumulates across cycles
- this makes the residual drift with cycle index
- the selected peak smears out
- `sigma` gets worse, selected fraction may drop, and merit falls

So the slope panel can look a bit like a differential version of the merit peak even though the code computes it independently.

In short:

- merit asks: how strong and narrow is the selected cluster?
- drift slope asks: does that selected cluster walk with cycle index?

Both are useful, but they answer different questions.

## Why Wings Appear Around The Main Peak Or Dip

When scanning trial periods, the X axis is the tested period value.

At each X point, the code folds all selected times with that trial period and asks how coherent the folded cluster is.

Near the true period, events line up strongly and produce the main peak or dip (depending on the metric).

If the trial period is slightly off, each next RF cycle is shifted a bit more than the previous one.
Across a finite number of cycles, this creates alternating partial alignment and partial cancellation.

That alternating pattern produces side lobes (the visible wings) around the central feature.
Because the alignment becomes weaker farther away, the side lobes decrease in size, which looks like a damped oscillation.

This expected finite-cycle sidelobe shape is often described by the Dirichlet kernel.

So in this scan, the wings are not a strange artifact by themselves.
They are the expected coherence pattern when you sweep trial period on the X axis and fold a finite cycle train.

## What Happens After The Best Period Is Chosen

1. Keep the best valid candidate, usually the one with the highest merit.
2. Build the final 1D best-phase profile.
3. Build the final `phase vs ToT` histogram.
4. Build cycle-residual summaries for the best-period solution.

Important:

- The search can use a subsample.
- The final `phase vs ToT` output is filled from all stored accepted pulses.

## How To Read The Best-Cycle Residual Plots

After the best period is chosen, the tool can also summarize the selected events by cycle index.

### Residual Versus Cycle Index

This plot uses:

- X axis: relative cycle index inside the `ch0`-referenced analysis window
- Y axis: centered residual of each selected event

So this is not a plot of absolute event time.
It is a plot of the event offset from its nearest expected best-period RF-cycle center.

The red mean line shows the average residual in each populated cycle.

### Local Slope Versus Cycle Index

This plot applies a moving local fit to the cycle-averaged residuals.

It shows where the residual trend is locally:

- rising with cycle index
- flat with cycle index
- falling with cycle index

This local slope plot is different from the single global slope number in the scan table.

- the global slope is one fitted number for the whole selected cluster
- the local slope plot shows how that trend changes across the window

## What This Method Assumes

This method assumes:

- one dominant folded phase cluster should still win after the full-phase seed-window search
- one global trial period can describe the data reasonably well
- one fixed peak window is enough to isolate the main cluster

So a good scan result does not automatically mean the final fold is artifact-free.
It only means that, under this model, one tested period performs best.

## Where To Look In The Code

- Data preparation and the search subsample: `finalize_window(...)`
- One trial-period evaluation: `evaluate_candidate(...)`
- Coarse-to-fine scan control: `run_multistage_scan(...)`
- Final best-period phase profile: `build_best_phase_profile(...)`
- Final best-period phase-vs-ToT output: `build_best_phase_tot_histogram(...)`