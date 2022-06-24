import os
import subprocess
import pandas as pd
import numpy as np
from itertools import combinations
from diffcomp.calder import CalderSubCompartments
from diffcomp.diff import CalderDifferentialCompartments, \
                          CalderRecursiveDifferentialSegmentator
from diffcomp.metrics import MeasureOfConcordance

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colorbar import ColorbarBase
import seaborn as sns


# ------------------------------------- #
#       BENCHMARKING OF DIFFCOMP        #
# ------------------------------------- #

# ********** #
# DATA SETUP #
# ********** #

# Data and Results directories
data_path = "data"
calder_path = os.path.join(data_path, "calder")
cores_path = os.path.join(data_path, "cores")
analysis_path = os.path.join(data_path, "analysis")

os.makedirs(calder_path, exist_ok=True)
os.makedirs(cores_path, exist_ok=True)
os.makedirs(analysis_path, exist_ok=True)

# Getting CALDER files from the WGD study
# We are going to focus on the 20 week post-WGD tumors to study the algorithm
# since they show higher amount of compartment repositioning
FIND_COMP_AND_CONCAT_CMD = "find {comp_path}/{sample}/CALDER_KR_50000 -name '*_sub_compartments.bed' | xargs cat > {outpath}"

compartments_path = "/mnt/etemp/luca/Projects/WorkInProgress/wgd/results/topology/compartments_CALDER"

samples = {
    # Mega maps
    "TP53_20wTumo1_30": "RPE_TP53_20w0T1_mega",
    "TP53_20wTumo2_30": "RPE_TP53_20w0T2_mega",
    "TP53_20wTumo3_30": "RPE_TP53_20w0T3_mega",
    "RPE_Control": "RPE_TP53_Ctrl_mega",

    # Replicas
    "RPE_C1_30": "RPE_TP53_Ctrl_1",
    "RPE_C2_30": "RPE_TP53_Ctrl_2",
    "TP53_20wTumo11_30": "RPE_TP53_20w0T1_1",
    "TP53_20wTumo12_30": "RPE_TP53_20w0T1_2",
    "TP53_20wTumo21_30": "RPE_TP53_20w0T2_1",
    "TP53_20wTumo22_30": "RPE_TP53_20w0T2_2",
    "TP53_20wTumo31_30": "RPE_TP53_20w0T3_1",
    "TP53_20wTumo32_30": "RPE_TP53_20w0T3_2",
}

print("----------------")
print("IMPORTING SAMPLE")
print("----------------")
for sample, sample_name in samples.items():
    print(sample)
    outpath = os.path.join(calder_path, f"{sample_name}_CompartmentDomains.bed")
    if not os.path.isfile(outpath):
        subprocess.run(FIND_COMP_AND_CONCAT_CMD.format(
                            comp_path = compartments_path,
                            sample = sample,
                            outpath = outpath),
                       shell=True)
        df = pd.read_csv(outpath, sep="\t", header=None)
        df[0] = df[0].map(lambda x: x if x != "chrNA" else "chrX")
        df.to_csv(outpath, sep="\t", header=False, index=False)


# ************************** #
# CORES ALGORITHM PARAMETERS #
# ************************** #

# Studying how CoRE change by changing the various parameters

BINSIZE=50000

GROUPS = {
    "RPE_TP53_Ctrl": ["RPE_TP53_Ctrl_mega", "RPE_TP53_Ctrl_1", "RPE_TP53_Ctrl_2"],
    "RPE_TP53_20w0T1": ["RPE_TP53_20w0T1_mega", "RPE_TP53_20w0T1_1", "RPE_TP53_20w0T1_2"],
    "RPE_TP53_20w0T2": ["RPE_TP53_20w0T2_mega", "RPE_TP53_20w0T2_1", "RPE_TP53_20w0T2_2"],
    "RPE_TP53_20w0T3": ["RPE_TP53_20w0T3_mega", "RPE_TP53_20w0T3_1", "RPE_TP53_20w0T3_2"],
}

SAMPLE_TO_GROUP = {
    v:k for k, x in GROUPS.items() for v in x
}

COMPARISONS = [
    ("RPE_TP53_Ctrl_mega", "RPE_TP53_20w0T1_mega")
]


# We study how CORES change by changing the maximum standard deviation of
# segments that is allowed
MIN_STDS = (np.arange(1, 10, dtype=int)/100).tolist() + \
           (np.arange(10, 20, dtype=int)/100).tolist() + \
           (np.arange(2, 10)/10).tolist()



print()
print("-"*110)
print("{:30}{:30}{:50}".format("SAMPLE1", "SAMPLE2", "STAGE"))
print("-"*110)
all_cores = {}
# We iterate for each sample comparison
for sample1, sample2 in COMPARISONS:

    # Building the set of control samples for this comparison
    # which is simply all the pairwise comparisons between samples
    # of the same group
    s1_replicas = GROUPS[SAMPLE_TO_GROUP[sample1]]
    s1_replicas = [CalderSubCompartments(path = os.path.join(calder_path, f"{x}_CompartmentDomains.bed")) for x in s1_replicas]
    s1_comparisons = list(combinations(s1_replicas, 2))
    s2_replicas = GROUPS[SAMPLE_TO_GROUP[sample2]]
    s2_replicas = [CalderSubCompartments(path = os.path.join(calder_path, f"{x}_CompartmentDomains.bed")) for x in s2_replicas]
    s2_comparisons = list(combinations(s2_replicas, 2))
    control_comparisons = s1_comparisons + s2_comparisons

    # We instantiate the segmentation algorithm
    segmentator = CalderRecursiveDifferentialSegmentator(binSize = BINSIZE)
    print(f"{sample1:30}{sample2:30}{'Building control distribution':50}")

    # We build the control distribution
    segmentator.build_control_distribution(control_comparisons)
    s1_comp = CalderSubCompartments(path = os.path.join(calder_path, f"{sample1}_CompartmentDomains.bed"))
    s2_comp = CalderSubCompartments(path = os.path.join(calder_path, f"{sample2}_CompartmentDomains.bed"))

    # We iterate on every value of sigma
    for min_std in reversed(MIN_STDS):
        minStd_str = f"Min. segment std = {min_std}"
        print(f"{sample1:30}{sample2:30}{minStd_str:50}")

        segments_path = os.path.join(cores_path, f"{sample1}_vs_{sample2}_CoREs_binSize_{BINSIZE}_minStd_{min_std}.tsv")
        if not os.path.isfile(segments_path):
            segmentator.minStd = min_std
            # Calling cores
            segments = segmentator.segment(s1_comp, s2_comp)
            segments.to_tsv(segments_path)
        else:
            segments = CalderDifferentialCompartments(path = segments_path)
        segments.to_bed(segments_path[:-4] + ".bed")
        all_cores[(sample1, sample2, BINSIZE, min_std)] = segments



# Studying the thresholds on p-value and fold-change

def do_pvalue_swipe(segments, pvalues):
    stat = []
    for p in pvalues:
        stat.append({
            "max_pval": p,
            "n_cores": segments.segmentation.assign(sign = lambda x: x.pvalue < p)["sign"].sum()
            })
    stat = pd.DataFrame.from_dict(stat)
    return stat

PVAL_SWIPE = [1, 0.5, 0.1, 0.05, 0.01, 2e-3,  1e-3, 5e-4, 2e-4, 1e-4]

all_segments_pval_swipe = []
for k, segments in all_cores.items():
    sample1, sample2, binsize, min_std = k
    segments_pval_swipe = do_pvalue_swipe(segments, PVAL_SWIPE)
    segments_pval_swipe['sample1'] = sample1
    segments_pval_swipe['sample2'] = sample2
    segments_pval_swipe['binsize'] = binsize
    segments_pval_swipe['min_std'] = min_std
    all_segments_pval_swipe.append(segments_pval_swipe)
all_segments_pval_swipe = pd.concat(all_segments_pval_swipe, axis=0, ignore_index=True)

for k, df in all_segments_pval_swipe.groupby(['sample1', 'sample2', 'binsize']):
    fig, ax = plt.subplots(1, 2, gridspec_kw={'width_ratios': [8, 1]})
    sns.lineplot(data = df[df.min_std <= 0.2],#.assign(min_std = lambda x: x.min_std.astype(str)),
                 x = "max_pval",
                 y = "n_cores",
                 hue = "min_std",
                 hue_norm = LogNorm(),
                 # hue_order = ["0.01"],
                 palette = "cool",
                 legend=False,
                 ax=ax[0])
    ax[0].set_box_aspect(1)
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("Max p-value")
    ax[0].set_ylabel("N. CoREs")
    sns.despine(ax = ax[0])
    cb = ColorbarBase(ax[1],
                      cmap = plt.cm.cool,
                      label = "Maximum standard deviation ($\sigma^*$)",
                      norm = LogNorm(vmin = df.min_std.min(), vmax = df.min_std.max()))
    fig.savefig(os.path.join(analysis_path, f"{k[0]}-{k[1]}-{k[2]}_pvalue_swipe.pdf"))
    plt.close(fig)


n_sign_cores_pval_swipe = []
for p in PVAL_SWIPE:
    n_sign_cores_pval_swipe.append(
        all_cores.assign(sign = lambda x: x.pvalue < p)\
            .groupby("comparison_name")['sign'].sum()\
            .to_frame("n_sign")\
            .reset_index()\
            .assign(p = p)
    )
n_sign_cores_pval_swipe = pd.concat(n_sign_cores_pval_swipe, axis=0, ignore_index=True)


MAX_PVAL = 0.01

fig = plt.figure()
sns.lineplot(data=n_sign_cores_pval_swipe, x = 'p', y = 'n_sign', hue = 'comparison_name')
plt.axvline(MAX_PVAL, color = 'black', linestyle = '--')
plt.xscale("log")
plt.yscale("symlog")
plt.ylim(0, 10000)
plt.xlabel("CoRE p-value ($p$)")
plt.ylabel("N. CoREs with p-value $\leq p$")
sns.despine(trim=False)
plt.legend(bbox_to_anchor=(1, 1), title = 'Comparison')
fig.savefig(os.path.join(results_path, "stats", "n_sign_cores_pval_swipe.pdf"))
plt.close(fig)



# Measuring similarities between CORES calls
measures_path = os.path.join(analysis_path, "measures.tsv")
if not os.path.isfile(measures_path):
    moc = MeasureOfConcordance()
    measures = []
    keys = list(all_cores.keys())

    print()
    print("-"*120)
    print("{:60}{:60}".format("CORES1", "CORES2"))
    print("-"*120)
    for i in range(len(keys)):
        ki = keys[i]
        for j in range(i, len(keys)):
            kj = keys[j]
            print("{:60}{:60}".format("-".join(map(str, ki)), "-".join(map(str, kj))))
            moc_ij = moc.compare(all_cores[ki].segmentation, all_cores[kj].segmentation, by_chrom=True)
            moc_ij = pd.Series(moc_ij).to_frame("moc").reset_index().rename(columns = {"index": "chr"})

            moc_ij['sample1_1'] = ki[0]
            moc_ij['sample2_1'] = ki[1]
            moc_ij['binsize_1'] = ki[2]
            moc_ij['min_std_1'] = ki[3]

            moc_ij['sample1_2'] = kj[0]
            moc_ij['sample2_2'] = kj[1]
            moc_ij['binsize_2'] = kj[2]
            moc_ij['min_std_2'] = kj[3]
            measures.append(moc_ij)
    measures = pd.concat(measures, axis=0, ignore_index=True)
    measures.to_csv(measures_path, sep="\t", index=False, header=True)
else:
    measures = pd.read_csv(measures_path, sep="\t")

for k, df in measures.groupby(["sample1_1", "sample2_1", "binsize_1", "sample1_2", "sample2_2", "binsize_2"]):
    kname = "-".join(map(str, k))
    H = pd.pivot_table(df, index = "min_std_1", columns = "min_std_2", values = "moc", aggfunc='mean')

    fig = plt.figure(figsize=(9, 9))
    sns.heatmap(H.sort_index(ascending=False),
                cmap='Reds',
                linewidth=0.5,
                xticklabels = False,
                cbar_kws = dict(location="top",
                                use_gridspec=False,
                                shrink=.5,
                                aspect=10,
                                label = "Measure of Concordance"),
                square=True)
    plt.yticks(rotation=0)
    plt.ylabel("Maximum standard deviation ($\sigma^*$)")
    plt.xlabel("")
    fig.savefig(os.path.join(analysis_path, f"{kname}_minStd_vs_MoC.pdf"))
    plt.close(fig)
