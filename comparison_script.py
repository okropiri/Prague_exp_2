import pandas as pd
import numpy as np

def analyze_tsv(file_path):
    # channel is at index 3, tot_ns is at index 8 (0-based)
    # Using pandas for potentially faster loading of large TSV
    df = pd.read_csv(file_path, sep='\t', comment='#')
    df_ch2 = df[df['channel'] == 2]
    
    total_count = len(df_ch2)
    max_tot = df_ch2['tot_ns'].max()
    tots = df_ch2['tot_ns'].values
    
    thresholds = [20, 30, 40, 60, 80, 100, 128]
    ge_counts = {f">={t}": np.sum(tots >= t) for t in thresholds}
    ge_fractions = {f">={t}": ge_counts[f">={t}"] / total_count if total_count > 0 else 0 for t in thresholds}
    
    bands = [(0, 20), (20, 40), (40, 60), (60, 80), (80, 100), (100, 128)]
    band_counts = {f"{low}-{high}": np.sum((tots >= low) & (tots < high)) for low, high in bands}
    band_counts[">=128"] = np.sum(tots >= 128)
    
    band_fractions = {k: v / total_count if total_count > 0 else 0 for k, v in band_counts.items()}
    
    return total_count, max_tot, ge_counts, ge_fractions, band_counts, band_fractions

def analyze_hist(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    metadata = {}
    data_start = 0
    for i, line in enumerate(lines):
        if line.startswith('#'):
            if '=' in line:
                key, val = line[1:].strip().split('=', 1)
                metadata[key] = val
        else:
            data_start = i
            break
            
    y_bins = int(metadata.get('y_bins', 0))
    y_max = float(metadata.get('y_max', 0))
    y_min = float(metadata.get('y_min', 0))
    bin_width = (y_max - y_min) / y_bins if y_bins > 0 else 0
    
    # x_bin_index y_bin_index count
    data = []
    for line in lines[data_start:]:
        parts = line.split()
        if len(parts) == 3:
            data.append([int(parts[0]), int(parts[1]), float(parts[2])])
    
    df = pd.DataFrame(data, columns=['x', 'y', 'count'])
    # y_bin_index is 1-based usually
    df['tot'] = (df['y'] - 0.5) * bin_width + y_min
    
    total_count = df['count'].sum()
    
    # We use the bin indices directly if they map well, but let's use the 'tot' calculated
    # Alternatively, use the bin index to be more precise if threshold aligns with bin borders
    # For these files, y_max=128, y_bins=512 or 256. 128/512 = 0.25 ns/bin. 128/256 = 0.5 ns/bin.
    # Thresholds 20, 30, 40, 60, 80, 100, 128 are exact multiples of 0.25 and 0.5.
    
    def get_sum(low, high=None):
        if high is None: # >= low
            mask = df['tot'] >= (low - 1e-7)
        else:
            mask = (df['tot'] >= (low - 1e-7)) & (df['tot'] < (high - 1e-7))
        return df.loc[mask, 'count'].sum()

    thresholds = [20, 30, 40, 60, 80, 100, 128]
    ge_counts = {f">={t}": get_sum(t) for t in thresholds}
    
    bands = [(0, 20), (20, 40), (40, 60), (60, 80), (80, 100), (100, 128)]
    band_counts = {f"{low}-{high}": get_sum(low, high) for low, high in bands}
    band_counts[">=128"] = get_sum(128)
    
    return total_count, y_max, ge_counts, band_counts

tsv_path = "/home/dachi/data6/Data/Refined_Data/NCAL_20us_Pos_3.4m_0000_rawRefined/NCAL_20us_Pos_3.4m_0000_rawRefined_pulses.tsv"
old_hist_path = "/home/dachi/data6/Dogma_analysis_by_Dachi/Results/NCAL_20us_Pos_3.4m_0000_ch0Ref_rfPhaseScan_ncalOnly_37to41/NCAL_20us_Pos_3.4m_0000_ch0Ref_rfPhaseScan_ncalOnly_37to41_best_phase_tot_hist.txt"
clean_hist_path = "/home/dachi/data6/Dogma_analysis_by_Dachi/Results/NCAL_20us_Pos_3.4m_0000_fullPhaseSeed/RF_period_scan/NCAL_20us_Pos_3.4m_0000_fullPhaseSeed_rf_period_scan_best_phase_tot_hist.txt"

print("--- TSV Analysis (Channel 2) ---")
tc, mtot, gec, gef, bc, bf = analyze_tsv(tsv_path)
print(f"Total Count: {tc}")
print(f"Max ToT: {mtot}")
for k in gec: print(f"{k}: {gec[k]} ({gef[k]:.6f})")
for k in bc: print(f"{k}: {bc[k]} ({bf[k]:.6f})")

print("\n--- Old Hist Analysis ---")
tc_o, ymax_o, gec_o, bc_o = analyze_hist(old_hist_path)
print(f"Total Count: {tc_o}")
print(f"Y Max: {ymax_o}")
for k in gec_o: print(f"{k}: {gec_o[k]} ({gec_o[k]/tc_o:.6f})")
for k in bc_o: print(f"{k}: {bc_o[k]} ({bc_o[k]/tc_o:.6f})")

print("\n--- Clean Hist Analysis ---")
tc_c, ymax_c, gec_c, bc_c = analyze_hist(clean_hist_path)
print(f"Total Count: {tc_c}")
print(f"Y Max: {ymax_c}")
for k in gec_c: print(f"{k}: {gec_c[k]} ({gec_c[k]/tc_c:.6f})")
for k in bc_c: print(f"{k}: {bc_c[k]} ({bc_c[k]/tc_c:.6f})")
