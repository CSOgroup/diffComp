%load_ext autoreload
%autoreload 2

import logging
import os
import sys
from diffcomp import CalderSubCompartments
from diffcomp.diff import CalderRecursiveDifferentialSegmentator, CalderDifferentialCompartments
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd
import numpy as np
# from pybedtools.bedtool import BedTool

# # from statsmodels.nonparametric.smoothers_lowess import lowess
# from abc import ABC, abstractmethod

from matplotlib import colormaps

from scipy.stats import poisson, gamma
# from multiprocessing import Pool

logformat = "[%(asctime)s] %(levelname)s:%(name)s: %(message)s"
logging.basicConfig(level=logging.DEBUG, stream=sys.stdout, format=logformat, datefmt="%Y-%m-%d %H:%M:%S")

BINSIZE=50000
CHROMS = [f'chr{x}' for x in range(1, 23)] + ['chrX']

GROUPS = {
    "RPE_TP53_Ctrl": ["RPE_TP53_Ctrl_mega", "RPE_TP53_Ctrl_1", "RPE_TP53_Ctrl_2"],
    "RPE_TP53_20w0T1": ["RPE_TP53_20w0T1_mega", "RPE_TP53_20w0T1_1", "RPE_TP53_20w0T1_2"]
}

SAMPLE_TO_GROUP = {
    v:k for k, x in GROUPS.items() for v in x
}


SAMPLE1 = "RPE_TP53_Ctrl_mega"
SAMPLE2 = "RPE_TP53_20w0T1_mega"


compartment_domains_path = "data/calder"

# Loading all compartment domains
comps = {}
for f in filter(lambda x: x.endswith("_CompartmentDomains.bed"), os.listdir(compartment_domains_path)):
	name = f.split("_CompartmentDomains.bed")[0]
	f_comp = CalderSubCompartments(os.path.join(compartment_domains_path, f))
	comps[name] = f_comp


# Building the set of control samples for this comparison
# which is simply all the pairwise comparisons between samples
# of the same group
s1_comp = comps[SAMPLE1]
s1_replicas = [comps[c] for c in GROUPS[SAMPLE_TO_GROUP[SAMPLE1]]]
s1_comparisons = list(combinations(s1_replicas, 2))
s2_replicas = [comps[c] for c in GROUPS[SAMPLE_TO_GROUP[SAMPLE2]]]
s2_comp = comps[SAMPLE2]
s2_comparisons = list(combinations(s2_replicas, 2))
control_comparisons = s1_comparisons + s2_comparisons


segmentator = CalderRecursiveDifferentialSegmentator(binSize = BINSIZE, statistical_test='gamma')
# segmentator.build_control_distribution([(s1, s2) for s1, s2 in control_comparisons])
segmentator.build_control_distribution(control_comparisons)
segments = segmentator.segment(comps[SAMPLE1], comps[SAMPLE2], chroms = CHROMS)

old_segments = pd.read_csv("tests/data/diff/TP53_20wTumo1_30_vs_RPE_Control_compartment_changes_recursive_segmentation.tsv", sep="\t")




# COMPARING CBS WITH DIFFCOMP
fig, ax = plt.subplots(2, 1, figsize=(20, 3))

X = segments.get_chromosomes("chrX")
cmap = colormaps['coolwarm']
v = X.signal["delta_rank"].values
vs  = np.convolve(v, np.ones(4)/4, mode='same')


ax[0].plot(X.signal.start.values, v, '.', color = '#D78521', markersize=2)
for t in X.segmentation[X.segmentation["gamma_pvalue"] <= 0.01].itertuples():
    ax[0].axvspan(t.start - BINSIZE//2, t.end - BINSIZE//2, alpha=0.2, color = cmap(np.clip((t.value*5 + 1)/2, 0, 1)))
ax[0].axhline(0, color='black', linewidth=1, linestyle='--')
ax[0].set_ylabel("$\Delta rank$")
ax[0].set_title("DiffComp algorithm (GAMMA p-value <= 0.01)", loc = "left", size=15)
ax[0].set_xlabel(f"Genomic coordinate (bp)")

ax[1].plot(X.signal.start.values, v, '.', color = '#D78521', markersize=2)

for t in old_segments[(old_segments.pvalue <= 0.01) & (old_segments.chr == "chrX")].itertuples():
    ax[1].axvspan(t.start - BINSIZE//2, t.end - BINSIZE//2, alpha=0.2, color = cmap(np.clip((t.value*5 + 1)/2, 0, 1)))
ax[1].axhline(0, color='black', linewidth=1, linestyle='--')
ax[1].set_ylabel("$\Delta rank$")
ax[1].set_title("DiffComp algorithm (OLD VERSION) (DELTA p-value <= 0.01)", loc = "left", size=15)
ax[1].set_xlabel(f"Genomic coordinate (bp)")


plt.show()



sns.histplot(data = segmentator.control_distribution['gamma']['chr10'], stat = 'density')
shape, loc, scale = segmentator.control_distribution['gamma-fit']['chr10']
x = np.linspace(gamma.ppf(0.001, shape, loc=loc, scale=scale), gamma.ppf(0.999, shape, loc=loc, scale=scale), 100)
plt.plot(x, gamma.pdf(x, shape, loc=loc, scale=scale), 'r-', label='gamma pdf')
plt.show()








segmentator_Delta = CalderRecursiveDifferentialSegmentator(binSize = BINSIZE, statistical_test="delta")
segmentator_Delta.build_control_distribution(control_comparisons)
segments_Delta = segmentator_Delta.segment(comps[SAMPLE1], comps[SAMPLE2])



sns.lmplot(data = segments.segmentation.assign(
						n_bins = lambda x: (x.end - x.start)//BINSIZE,
						abs_value = lambda x: x.value.abs(),
						negLog10Pvalue = lambda x: -1*np.log10(x["gamma_pvalue"])
					), x = "n_bins", y = "gamma_statistic", hue = 'chr', palette = "coolwarm_r")
plt.legend().set_visible(False)
plt.show()


sns.scatterplot(data=segments.segmentation.assign(
						n_bins = lambda x: (x.end - x.start)//BINSIZE,
						negLog10Pvalue = lambda x: -1*np.log10(x["gamma_pvalue"])
					),
				x = "n_bins", y = "negLog10Pvalue", size = "n_bins", hue = 'chr')
plt.show()



sns.displot(data = segments.segmentation.assign(
						n_bins = lambda x: (x.end - x.start)//BINSIZE,
						negLog10Pvalue = lambda x: -1*np.log10(x["gamma_pvalue"])
					), x = "negLog10Pvalue", hue = "chr", kind='ecdf', palette="coolwarm_r")
plt.show()



sns.displot(data = segments.segmentation.assign(
						n_bins = lambda x: (x.end - x.start)//BINSIZE,
						negLog10Pvalue = lambda x: -1*np.log10(x.pvalue)
					),)



sns.scatterplot(data=segments_Delta.segmentation.assign(
						n_bins = lambda x: (x.end - x.start)//BINSIZE,
						negLog10Pvalue = lambda x: -1*np.log10(x.pvalue)
					),
				x = "value", y = "negLog10Pvalue", size = "n_bins")
plt.show()




X = pd.concat([
	segmentator._null_dist.assign(data = 'control')[['chr', "start", 'end', "abs_value", 'n_bins', 'length', 'data']],
	segments.segmentation.assign(length = lambda x: x.end - x.start, abs_value = lambda x: x.value.abs(), n_bins = lambda x: x.length // BINSIZE, data = "real")\
				[['chr', "start", 'end', "abs_value", 'n_bins', 'length', "data"]]
], axis=0, ignore_index=True)


X['perc_bins'] = X['n_bins'] / (segments.segmentation.end - segments.segmentation.start).sum()

sns.displot(data = X, x = "n_bins", hue = "data", kind='kde', common_norm=False)
plt.show()



# Generate some sample data
data = gamma.rvs(a=5, size=1000)

# Fit the data to a gamma distribution
shape, loc, scale = gamma.fit(data)


# null_dist = segments.
# Create a histogram of the data
plt.hist(data, bins=50, alpha=0.5, density=True)

# Create a range of x-values for the fitted gamma distribution
x = np.linspace(gamma.ppf(0.001, shape, loc=loc, scale=scale), gamma.ppf(0.999, shape, loc=loc, scale=scale), 100)

# Plot the fitted gamma distribution
plt.plot(x, gamma.pdf(x, shape, loc=loc, scale=scale), 'r-', label='gamma pdf')

# Add a legend to the plot
plt.legend()

# Show the plot
plt.show()




plt.hist(null_dist, bins=30, density=True)
plt.show()







X['formula'] = (X['n_bins']*X['abs_value'])

null_dist = X.loc[X.data == 'control', "formula"].values

fit_alpha, fit_loc, fit_beta=gamma.fit(null_dist)





X['pvalue'] = X.formula.map(lambda x: ((null_dist > x).sum() + 1) / (null_dist.shape[0] + 1) )
X['negLog10Pvalue'] = -1*np.log10(X.pvalue)

sns.displot(data = X, x = 'formula', hue = "data", kind='hist', stat = "probability", common_norm=False, hue_order=['control'])
x = np.arange(0, 5, 0.01)
y = gamma.pdf(x, fit_alpha, fit_loc, fit_beta)
plt.plot(x, y)
plt.show()

sns.displot(data = X, x = 'formula', hue = "data", kind='ecdf', hue_order=['control'])
plt.show()





seg_n_bins = 10
(segmentator._null_dist.n_bins - seg_n_bins).abs().sort_values().index[:]




segments = segmentator.segment(control_comparisons[0][0].get_chromosomes('chr19'), control_comparisons[0][1].get_chromosomes('chr19'))



