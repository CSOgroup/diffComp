from importlib import reload
import diffcomp
import os
import pandas as pd
import numpy as np
from pybedtools.bedtool import BedTool
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import seaborn as sns


def get_s12_rank(comps1, comps2, binSize=50000):
	segmentator = diffcomp.CalderRecursiveDifferentialSegmentator(binSize, 0.1)
	s12_rank = segmentator.get_binned_deltaRank(comps1, comps2)
	return s12_rank



def closest(bounds1, bounds2):
	"""Find the closest boundary in bounds2 to the boundaries in bounds1.
	The two lists of boundaries have to be pre-sorted"""
	if isinstance(bounds1, list):
		bounds1 = np.array(bounds1)
	if isinstance(bounds2, list):
		bounds2 = np.array(bounds2)
	closest = bounds2[np.abs(bounds1.reshape(-1, 1) - bounds2.reshape(1, -1)).argmin(axis=1)]
	return closest

def get_valid_bins(domains, binsize):
	"""Gets all bins covered by a domain"""
	valid_positions = []
	for d in domains.itertuples():
		dv = np.arange(d.start, d.end + binsize, binsize)
		if (len(valid_positions) > 0):
			if valid_positions[-1][-1] == dv[0]:
				dv = dv[1:]
		valid_positions.append(dv)
	valid_positions = np.concatenate(valid_positions)
	return valid_positions

compartment_domains_path = "data/calder"
output_path = "data/analysis/new_implementation"
os.makedirs(output_path, exist_ok=True)

def get_sign_change_positions(values):
	# Get the signs of values. Positive = 1, Negative = -1. NaN values will be assigned 0.
	signs = np.nan_to_num(np.sign(values))

	# Get positions of sign change
	signchange = ((np.roll(signs, 1) - signs) != 0).astype(int)
	signchange[0] = 0
	return signchange

# Loading all compartment domains
comps = {}
for f in filter(lambda x: x.endswith("_CompartmentDomains.bed"), os.listdir(compartment_domains_path)):
	name = f.split("_CompartmentDomains.bed")[0]
	f_comp = diffcomp.CalderSubCompartments(os.path.join(compartment_domains_path, f))
	comps[name] = f_comp


SAMPLES = [
	"RPE_TP53_Ctrl_1",
	"RPE_TP53_Ctrl_2",
	"RPE_TP53_Ctrl_mega",
	"RPE_TP53_20w0T1_mega",
	"CPA_TP53_Ctrl_c3_mega"
]

COMPARISONS = {
	"rep1 vs rep2"			: ("RPE_TP53_Ctrl_1", "RPE_TP53_Ctrl_2"),
	"control vs tumor"		: ("RPE_TP53_Ctrl_mega", "RPE_TP53_20w0T1_mega"),
	"RPE vs CP-A"			: ("RPE_TP53_Ctrl_mega", "CPA_TP53_Ctrl_c3_mega")
}

BINSIZE = 50000
CHROM = "chr5"


EXAMPLE_REGION = (10000000, 30000000)

## Plotting an example of Calder segmentation in the following samples:
## - Replicate 1
## - Replicate 2
## - Control
## - Tumor
## - Different cell line

fig, ax = plt.subplots(len(SAMPLES), 1,
					   figsize = (20, len(SAMPLES)*2),
					   tight_layout=True)
for i, sample in enumerate(SAMPLES):
	domains = comps[sample].get_chromosome(CHROM).domains.sort_values(['chr', 'start', 'end'])

	starts = domains.start.values
	ends = domains.end.values
	rank = domains.domain_rank.values

	cmap = get_cmap("coolwarm")
	colors = cmap(rank)

	ax[i].hlines(rank,
				 starts,
				 ends,
				 colors = colors,
				 linewidth=3)
	ax[i].set_xlim(*EXAMPLE_REGION)
	ax[i].set_ylim(0, 1)
	ax[i].set_title(sample, fontsize=20)
	ax[i].set_ylabel("Rank")
	xticks = ax[i].get_xticks()
	ax[i].set_xticks(xticks, labels = [f"{int(x/1e6)} Mb" for x in xticks])
	ax[i].spines['top'].set_visible(False)
	ax[i].spines['right'].set_visible(False)
fig.savefig(os.path.join(output_path, "sample_ranks.pdf"))
plt.close(fig)


## Show the derivative of ranks for each sample
fig, ax = plt.subplots(len(SAMPLES), 1,
					   figsize = (20, len(SAMPLES)*2),
					   tight_layout=True)
for i, sample in enumerate(SAMPLES):
	bins = comps[sample].get_chromosome(CHROM).binnify(BINSIZE).sort_values(['chr', 'start', 'end'])

	starts = bins.start.values
	ends = bins.end.values
	rank = bins.domain_rank.values
	der_rank = np.diff(rank, append=rank[-1])

	ax[i].hlines(der_rank,
				 starts,
				 ends,
				 colors = "black",
				 linewidth=3)
	ax[i].set_xlim(*EXAMPLE_REGION)
	ax[i].set_ylim(-1, 1)
	ax[i].set_title(sample, fontsize=20)
	ax[i].set_ylabel("$\partial$Rank")
	xticks = ax[i].get_xticks()
	ax[i].set_xticks(xticks, labels = [f"{int(x/1e6)} Mb" for x in xticks])
	ax[i].spines['top'].set_visible(False)
	ax[i].spines['right'].set_visible(False)
fig.savefig(os.path.join(output_path, "sample_der_ranks.pdf"))
plt.close(fig)



## Extract Delta-Rank vectors for each comparison
all_s12_rank = []
for comparison, samples in COMPARISONS.items():
	s12_rank = get_s12_rank(comps[samples[0]].get_chromosome(CHROM),
							comps[samples[1]].get_chromosome(CHROM))\
				.assign(comparison = comparison)
	all_s12_rank.append(s12_rank)
all_s12_rank = pd.concat(all_s12_rank, axis=0, ignore_index=True)


# Show differences in Delta-Rank between different comparison types
fig, ax = plt.subplots(len(COMPARISONS), 1,
					   figsize = (20, len(COMPARISONS)*2),
					   tight_layout=True)
i = 0
for comparison in COMPARISONS.keys():
	df = all_s12_rank[all_s12_rank.comparison == comparison]
	starts = df.start.values
	ends = df.end.values
	values = df.delta_rank.values
	# sign_changes = get_sign_change_positions(values)

	ax[i].hlines(values,
				 starts,
				 ends,
				 linewidth=3)
	ax[i].axhline(0, color='black', linewidth=1, linestyle='--')
	# for p in np.where(sign_changes > 0)[0]:
	# 	ax[i].axvline(p, color = 'red', linewidth=0.5, linestyle='--')
	ax[i].set_xlim(*EXAMPLE_REGION)
	ax[i].set_title(comparison, fontsize=20)
	ax[i].set_ylabel("$\Delta$rank")
	xticks = ax[i].get_xticks()
	ax[i].set_xticks(xticks, labels = [f"{int(x/1e6)} Mb" for x in xticks])
	ax[i].set_ylim(-1, 1)
	ax[i].text(1, 1.1, CHROM, transform=ax[i].transAxes, ha='right', fontsize=15)
	i += 1
fig.savefig(os.path.join(output_path, "comparisons_deltaRank.pdf"))
plt.close(fig)


# Show differences between derivatives of rank for each comparison
all_s12_der_rank = []
fig, ax = plt.subplots(len(COMPARISONS), 1,
					   figsize = (20, len(COMPARISONS)*2),
					   tight_layout=True)
i = 0
for comparison, samples in COMPARISONS.items():
	bins1 = comps[samples[0]].get_chromosome(CHROM).binnify(BINSIZE).sort_values(['chr', 'start', 'end'])
	bins2 = comps[samples[1]].get_chromosome(CHROM).binnify(BINSIZE).sort_values(['chr', 'start', 'end'])
	bins = bins1[['chr', 'start', 'end', 'domain_rank']]\
				.rename(columns = {'domain_rank': 'domain_rank_1'})\
				.merge(bins2[['chr', 'start', 'end', 'domain_rank']]\
							.rename(columns = {'domain_rank': 'domain_rank_2'}),
					   on = ['chr', 'start', 'end'])\
				.sort_values(['chr', 'start', 'end'])\
				.reset_index(drop=True)
	bins['delta_rank'] = bins.domain_rank_2 - bins.domain_rank_1
	bins['der_domain_rank_1'] = np.diff(bins.domain_rank_1.values, append=bins.domain_rank_1.values[-1])
	bins['der_domain_rank_2'] = np.diff(bins.domain_rank_2.values, append=bins.domain_rank_2.values[-1])
	bins['delta_der_domain_rank'] = bins['der_domain_rank_2'] - bins['der_domain_rank_1']
	bins['comparison'] = comparison
	all_s12_der_rank.append(bins)

	starts = bins.start.values
	ends = bins.end.values
	delta_rank = bins.delta_rank

	ax[i].hlines(delta_rank,
			 starts,
			 ends,
			 colors = 'green',
			 linewidth=3)

	nonZeroBins = bins[(bins.delta_der_domain_rank != 0) & (~bins.delta_der_domain_rank.isnull())]

	delta_der_rank = nonZeroBins.delta_der_domain_rank.values
	ax[i].hlines(delta_der_rank,
				 starts,
				 ends,
				 colors = 'grey',
				 linewidth=3)
	ax[i].axhline(0, color='black', linewidth=1, linestyle='--')
	# for p in np.where(sign_changes > 0)[0]:
	# 	ax[i].axvline(p, color = 'red', linewidth=0.5, linestyle='--')
	ax[i].set_xlim(*EXAMPLE_REGION)
	ax[i].set_title(comparison, fontsize=20)
	ax[i].set_ylabel("$\Delta\partial$Rank")
	xticks = ax[i].get_xticks()
	ax[i].set_xticks(xticks, labels = [f"{int(x/1e6)} Mb" for x in xticks])
	ax[i].set_ylim(-1, 1)
	ax[i].text(1, 1.1, CHROM, transform=ax[i].transAxes, ha='right', fontsize=15)
	i += 1
fig.savefig(os.path.join(output_path, "comparisons_deltaDerRank.pdf"))
plt.close(fig)
all_s12_der_rank = pd.concat(all_s12_der_rank, axis=0, ignore_index=True)






# Show differences in compartment domain boundaries between different comparison types
fig, ax = plt.subplots(len(SAMPLES), 1,
					   figsize = (20, len(SAMPLES)*2),
					   tight_layout=True)
for i, sample in enumerate(SAMPLES):
	domains = comps[sample].get_chromosome(CHROM).domains.sort_values(['chr', 'start', 'end'])
	starts = domains.start.values
	ends = domains.end.values
	ax[i].hlines(np.repeat(0, np.arange(0, starts.shape[0], 2).shape[0]),
			   starts[np.arange(0, starts.shape[0], 2)],
			   ends[np.arange(0, ends.shape[0], 2)])

	ax[i].hlines(np.repeat(1, np.arange(1, starts.shape[0], 2).shape[0]),
		   starts[np.arange(1, starts.shape[0], 2)],
		   ends[np.arange(1, ends.shape[0], 2)])
	ax[i].set_ylim(-2, 3)
	ax[i].set_xlim(10000000, 20000000)
	ax[i].set_title(sample, fontsize=20)
	ax[i].set_yticks([])
	xticks = ax[i].get_xticks()
	ax[i].set_xticks(xticks, labels = [f"{int(x/1e6)} Mb" for x in xticks])
	ax[i].spines['top'].set_visible(False)
	ax[i].spines['right'].set_visible(False)
	ax[i].spines['left'].set_visible(False)
fig.savefig(os.path.join(output_path, "samples_bounds.pdf"))
plt.close(fig)







sns.relplot(data = pd.concat([
		s12_rank.assign(name = "tumor vs control"),
		control12_rank.assign(name = "replicate vs replicate")
	], axis=0, ignore_index=True),
	x = "domain_rank_1",
	y = "delta_rank",
	col = "name")
plt.show()

control_bounds = comps_control.get_domain_boundaries(CHROM)
tumor_bounds = comps_tumor.get_domain_boundaries(CHROM)
rep1_bounds = comps_rep1.get_domain_boundaries(CHROM)
rep2_bounds = comps_rep2.get_domain_boundaries(CHROM)


fig = plt.figure()
plt.hlines(comps_control.domains.domain_rank, xmin=comps_control.domains.start.values, xmax = comps_control.domains.end.values)
for b in control_bounds:
	plt.axvline(b, linewidth=0.5, linestyle='--', color='black')
plt.xlim(0, 20000000)
plt.show()



def get_aligned_boundaries(bounds1, bounds2, binsize):
	close_bounds = closest(bounds1, bounds2)
	return pd.DataFrame({
			"boundary": bounds1,
			"closest_boundary": close_bounds,
			'distance': np.abs(bounds1 - close_bounds),
			'bin_distance': np.abs(bounds1 - close_bounds)//binsize
		})



aligned_bounds = get_aligned_boundaries(control_bounds, tumor_bounds, BINSIZE).assign(name = "tumor vs control")
control_aligned_bounds = get_aligned_boundaries(rep1_bounds, rep2_bounds, BINSIZE).assign(name = "replicate vs replicate")



random_bounds = np.sort(np.random.choice(np.arange(np.min(control_bounds),
												   np.max(control_bounds),
												   BINSIZE),
					 		size = control_bounds.shape[0],
					 		replace=False))
random_aligned_bounds = get_aligned_boundaries(control_bounds, random_bounds, BINSIZE).assign(name = "control vs random")

theoretical_random_mean_bin_distance = ((np.max(control_bounds) - np.min(control_bounds)) / control_bounds.shape[0])//(2*BINSIZE)

sns.displot(data = pd.concat([
		aligned_bounds,
		control_aligned_bounds,
		random_aligned_bounds
	], axis=0, ignore_index=True),
	x = "bin_distance",
	hue = "name",
	discrete=True)
plt.show()





sample1 = "HTB9_DMSO_utx-stag2"
sample2 = "HTB9-UTXko-Clone3_DMSO_utx-stag2"
compartment_domains_path = "data/calder"

sample1_comp = diffcomp.CalderSubCompartments(os.path.join(compartment_domains_path, f"{sample1}_CompartmentDomains.bed"),
											  genome = "hg19",
											  coordinates='one-based')

sample2_comp = diffcomp.CalderSubCompartments(os.path.join(compartment_domains_path, f"{sample2}_CompartmentDomains.bed"),
											  genome = "hg19",
											  coordinates='one-based')

bins = BedTool().window_maker(w = 100000, genome = 'hg19')
sample1_binned = bins.map(BedTool.from_dataframe(sample1_comp.domains).sort(),
						  c=[5, 6], o = ['distinct'], null = np.NaN)\
					 .to_dataframe(names = ['chr', 'start', 'end', 'compartment1', 'domainRank1'])
sample2_binned = bins.map(BedTool.from_dataframe(sample2_comp.domains).sort(),
						  c=[5, 6], o = ['distinct'], null = np.NaN)\
					 .to_dataframe(names = ['chr', 'start', 'end', 'compartment2', 'domainRank2'])

aligned_comp = sample1_binned.merge(sample2_binned, on =['chr', 'start', 'end'])
aligned_comp['deltaRank'] = aligned_comp.domainRank2 - aligned_comp.domainRank1
aligned_comp = aligned_comp.dropna()
