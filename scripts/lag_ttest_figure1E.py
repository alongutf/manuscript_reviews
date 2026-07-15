"""
Two-sample t-test on scanlag lag-time distributions (Figure 1E).

Unit of replication = biological replicate (n = 3 per condition), so we first
collapse each replicate to a single number (its mean lag time) and then run the
t-test on those three means. This is more conservative and more defensible than
pooling all colonies, because colonies within one plate are not independent.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
# Run from scripts/ so that the repo root is the parent of the working dir.
REPO = os.path.dirname(os.getcwd())
data_dir = os.path.join(REPO, "scanlag_data", "exp2")

# --------------------------------------------------------------------------
# Normalization (identical to what figure1.py panel_E plots)
# --------------------------------------------------------------------------
# Each colony's raw X value is its absolute appearance time. The plot expresses
# lag *relative* to the earliest exponential colony, so it subtracts a single
# constant t0 = min appearance time of the REP3 exponential culture.
# NOTE: because t0 is one constant subtracted from every group equally, it shifts
# all means by the same amount and therefore does NOT change any t-test result.
# We apply it only to stay faithful to the plotted data.
t0 = np.min(pd.read_csv(os.path.join(data_dir, "REP3EXP_t00Min_ax1.csv"))["X"])

# Map each condition label to the substring that identifies it in the filenames.
# Files look like: REP1SHX_t31756Min_ax1.csv, REP3CASP_..., REP4EXP_..., etc.
groups = {
    "Exponential": "EXP",    # exponential-phase control
    "Reg-Arrest":  "CASP",   # regulated arrest  (chloramphenicol)
    "Dis-Arrest":  "SHX",    # dysregulated arrest (serine hydroxamate)
}

# rep_means[condition][replicate] = mean lag time of that single plate
rep_means = {k: {} for k in groups}

# --------------------------------------------------------------------------
# Reconstruct per-colony lag times from each survival curve, then take the mean
# --------------------------------------------------------------------------
# Each CSV is an empirical survival function: columns X (appearance time) and Y
# (fraction of colonies not yet appeared). Y starts at 1 and drops by 1/N at each
# colony, so every row *after the first* corresponds to exactly one colony whose
# lag time is that row's X. Dropping the first row removes the Y = 1 baseline.
for f in sorted(os.listdir(data_dir)):
    d = pd.read_csv(os.path.join(data_dir, f))
    lag_times = d["X"].values[1:] - t0          # per-colony lag times (normalized)

    # Replicate id = the filename prefix before the condition tag (REP1/REP3/REP4)
    rep = f.split("EXP")[0].split("CASP")[0].split("SHX")[0]

    for condition, tag in groups.items():
        if tag in f:
            rep_means[condition][rep] = lag_times.mean()   # collapse plate -> 1 value

# --------------------------------------------------------------------------
# Report the per-replicate means (sanity check)
# --------------------------------------------------------------------------
print("Per-replicate mean lag time (min):")
for condition in groups:
    for rep, m in sorted(rep_means[condition].items()):
        print(f"  {condition:12s} {rep}: {m:7.1f}")
    print()

# --------------------------------------------------------------------------
# The actual test: 3 means vs 3 means
# --------------------------------------------------------------------------
def rep_ttest(a_dict, b_dict, name_a, name_b):
    a = np.array(list(a_dict.values()))   # 3 replicate means, group A
    b = np.array(list(b_dict.values()))   # 3 replicate means, group B

    # Welch's t-test (equal_var=False) does NOT assume equal variance -> safer default.
    # Student's t-test (equal_var=True) pools the variances; shown for comparison.
    tw, pw = ttest_ind(a, b, equal_var=False)
    ts, ps = ttest_ind(a, b, equal_var=True)

    print(f"{name_a} vs {name_b}")
    print(f"  {name_a}: mean={a.mean():.1f}, sd={a.std(ddof=1):.1f}, n={a.size}")
    print(f"  {name_b}: mean={b.mean():.1f}, sd={b.std(ddof=1):.1f}, n={b.size}")
    print(f"  Welch's   t = {tw:.3f}, p = {pw:.4f}")
    print(f"  Student's t = {ts:.3f}, p = {ps:.4f}  (df={a.size + b.size - 2})")
    print()

rep_ttest(rep_means["Dis-Arrest"], rep_means["Reg-Arrest"],  "Dis-Arrest", "Reg-Arrest")
rep_ttest(rep_means["Dis-Arrest"], rep_means["Exponential"], "Dis-Arrest", "Exponential")
