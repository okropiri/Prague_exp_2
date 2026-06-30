import numpy as np
import uproot
from uproot.writing.identify import to_TAxis, to_TH2x

def make_uniform_histogram2d(
    counts: np.ndarray,
    title: str,
    x_title: str,
    y_title: str,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
):
    x_bins = counts.shape[1]
    y_bins = counts.shape[0]
    storage = np.zeros((y_bins + 2, x_bins + 2), dtype=np.float64)
    storage[1:-1, 1:-1] = counts

    x_centers = x_min + (np.arange(x_bins, dtype=np.float64) + 0.5) * ((x_max - x_min) / x_bins)
    y_centers = y_min + (np.arange(y_bins, dtype=np.float64) + 0.5) * ((y_max - y_min) / y_bins)
    x_centers_2d = x_centers[np.newaxis, :]
    y_centers_2d = y_centers[:, np.newaxis]
    entries = float(counts.sum(dtype=np.float64))

    return to_TH2x(
        fName=None,
        fTitle=title,
        data=storage.ravel(),
        fEntries=entries,
        fTsumw=entries,
        fTsumw2=float(np.square(counts, dtype=np.float64).sum(dtype=np.float64)),
        fTsumwx=float((counts * x_centers_2d).sum(dtype=np.float64)),
        fTsumwx2=float((counts * np.square(x_centers_2d)).sum(dtype=np.float64)),
        fTsumwy=float((counts * y_centers_2d).sum(dtype=np.float64)),
        fTsumwy2=float((counts * np.square(y_centers_2d)).sum(dtype=np.float64)),
        fTsumwxy=float((counts * x_centers_2d * y_centers_2d).sum(dtype=np.float64)),
        fSumw2=None,
        fXaxis=to_TAxis("xaxis", x_title, x_bins, float(x_min), float(x_max)),
        fYaxis=to_TAxis("yaxis", y_title, y_bins, float(y_min), float(y_max)),
    )

# Test setup
x_bins, y_bins = 10, 5
counts = np.zeros((y_bins, x_bins))
# Use distinctive values at specific bins
counts[1, 2] = 42.0  # y_bin=1 (index 1), x_bin=2 (index 2)
counts[3, 7] = 17.0  # y_bin=3 (index 3), x_bin=7 (index 7)

filename = "test_root.root"
with uproot.recreate(filename) as f:
    f["h2"] = make_uniform_histogram2d(counts, "Test", "X", "Y", 0, 10, 0, 5)

# Verification
with uproot.open(filename) as f:
    h = f["h2"]
    values = h.values()
    # uproot.TH2.values() returns a 2D array where rows are X and columns are Y by default, 
    # OR it matches the storage. Actually it returns [x_bin, y_bin] or [y_bin, x_bin]?
    # Let's check the shape and content.
    print(f"Read shape: {values.shape}")
    print(f"Value at [2, 1]: {values[2, 1]}")
    print(f"Value at [7, 3]: {values[7, 3]}")
    
    if values[2, 1] == 42.0 and values[7, 3] == 17.0:
        print("SUCCESS: Indices match [x, y] access in uproot")
    else:
        print("FAILURE: Bin mismatch")
