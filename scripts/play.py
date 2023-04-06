from diffcomp.diff import CalderDifferentialSegmentator, CalderDifferentialCompartments
from diffcomp import CalderSubCompartments
import os
from typing import Generic, Optional, List, Tuple
from itertools import combinations
import pandas as pd
import numpy as np
from pybedtools.bedtool import BedTool
import logging
from statsmodels.nonparametric.smoothers_lowess import lowess
from abc import ABC, abstractmethod
import matplotlib.pyplot as plt
from matplotlib import colormaps
import seaborn as sns
from scipy.stats import poisson, gamma
from multiprocessing import Pool


BINSIZE=50000

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



class CalderDifferentialCompartments:
	CALDER_DIFFERENTIAL_COMPARTMENTS_HEADER = ['chr', 'start', 'end', 'value', "rank1", "comp1", "rank2", "comp2", 'statistic', 'pvalue']

	def __init__(self,
				 path: str = None,
				 input_df: Optional[pd.DataFrame] = None,
				 signal_path: Optional[str] = None,
				 signal_df: Optional[pd.DataFrame] = None,
				 expected_std: Optional[pd.DataFrame] = None):
		if (path is None) and (input_df is not None):
			self._check_header(input_df)
			self._segs = input_df
		elif (path is not None) and (input_df is None):
			self._segs = self.read(path)
		else:
			raise ValueError("'path' and 'input_df' are mutually exclusive, and at least one has to be provided")

		if (signal_path is None) and (signal_df is not None):
			self._signal = signal_df
		elif (signal_path is not None) and (signal_df is None):
			self._signal = self.read_signal(signal_path)
		elif (signal_path is None) and (signal_df is None):
			self._signal = None
		else:
			raise ValueError("'signal_path' and 'signal_df' are mutually exclusive")

		if expected_std is not None:
			self._expected_std = expected_std
		else:
			self._expected_std = None


	@classmethod
	def read(cls, path: str):
		segs = pd.read_csv(path, sep='\t')
		cls._check_header(segs)
		return segs

	@classmethod
	def read_signal(cls, path: str):
		signal = pd.read_csv(path, sep="\t")
		return signal

	@classmethod
	def _check_header(cls, segs: pd.DataFrame):
		if segs.columns[:len(cls.CALDER_DIFFERENTIAL_COMPARTMENTS_HEADER)].tolist() != cls.CALDER_DIFFERENTIAL_COMPARTMENTS_HEADER:
			raise ValueError("Calder differential segmentation file must have as initial columns " + ", ".join(cls.CALDER_DIFFERENTIAL_COMPARTMENTS_HEADER))

	@property
	def segmentation(self):
		return self._segs

	@property
	def signal(self):
		return self._signal

	@property
	def expected_std_by_segment_size(self):
		return self._expected_std

	@staticmethod
	def concat(segmentations: list[CalderDifferentialCompartments]) -> CalderDifferentialCompartments:
		all_segs = [x.segmentation for x in segmentations]
		all_segs = pd.concat(all_segs, axis=0, ignore_index=True)
		all_segs = all_segs.sort_values(['chr', 'start', 'end']).reset_index(drop=True)
		all_signals = [x.signal for x in segmentations if x.signal is not None]
		if len(all_signals) > 0:
			all_signals = pd.concat(all_signals, axis=0, ignore_index=True)
			all_signals = all_signals.sort_values(['chr', 'start', 'end']).reset_index(drop=True)
		else:
			all_signals = None
		all_expected_std = [x.expected_std_by_segment_size for x in segmentations if x.expected_std_by_segment_size is not None]
		if len(all_expected_std) > 0:
			all_expected_std = pd.concat(all_expected_std, axis=0, ignore_index=True)
		else:
			all_expected_std = None
		return CalderDifferentialCompartments(input_df=all_segs, signal_df = all_signals, expected_std=all_expected_std)

	def to_tsv(self, path: str):
		self._segs.to_csv(path, sep="\t", index=False, header=True)

	def to_bed(self, path: str):
		cmap = LinearSegmentedColormap.from_list("CoREs_colormap", ['blue', "lightgrey", "red"])
		self._segs.assign(
					  name = lambda x: x.chr + ":" + x.start.astype(str) + "-" + x.end.astype(str),
					  norm_value = lambda x: (x.value + 1)/2,
					  score = lambda x: x.value,
					  strand = ".",
					  thickStart = lambda x: x.start,
					  thickEnd = lambda x: x.end,
					  color = lambda x: list(map(lambda x: ",".join(x), (cmap(x.norm_value.values)*255).astype(int).astype(str)[:, :-1].tolist())))\
				[['chr', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'color']]\
				.to_csv(path, sep='\t', index=False, header=False)

	def write_signal(self, path: str):
		if self._signal is not None:
			self._signal.dropna().to_csv(path, sep="\t", index=False, header=False)


class CalderDifferentialSegmentator(ABC):
	@abstractmethod
	def segment_chromosome(self,
						   s1: CalderSubCompartments,
						   s2: CalderSubCompartments) -> CalderDifferentialCompartments:
		pass

	@property
	@abstractmethod
	def binSize(self) -> int:
		pass

	@property
	@abstractmethod
	def n_cpus(self) -> int:
		pass

	def get_binned_deltaRank(self, s1: CalderSubCompartments, s2: CalderSubCompartments) -> pd.DataFrame:
		s1_binned, s2_binned = s1.binnify(self._binSize), s2.binnify(self._binSize)
		s1_rank = s1_binned[['chr', 'start', 'end', 'domain_rank', 'compartment_label_8']]
		s2_rank = s2_binned[['chr', 'start', 'end', 'domain_rank', 'compartment_label_8']]
		s12_rank = s1_rank.merge(s2_rank, on=['chr', 'start', 'end'], suffixes=("_1", "_2"))
		s12_rank['delta_rank'] = s12_rank.domain_rank_2 - s12_rank.domain_rank_1
		return s12_rank

	# def segment(self,
	# 			s1: CalderSubCompartments,
	# 			s2: CalderSubCompartments) -> CalderDifferentialCompartments:
	# 	_logger = logging.getLogger(self.__class__.__name__)
	# 	_logger.info("Binning sample 1")
	# 	s1.binnify(self.binSize)
	# 	_logger.info("Binning sample 2")
	# 	s2.binnify(self.binSize)
	# 	chroms = set(s1.domains.chr).intersection(set(s2.domains.chr))
	# 	all_chrom_seg = []
	# 	for chrom in chroms:
	# 		_logger.info(f"Processing {chrom}")
	# 		chrom_seg = self.segment_chromosome(s1.get_chromosomes(chrom), s2.get_chromosomes(chrom))
	# 		all_chrom_seg.append(chrom_seg)
	# 	return CalderDifferentialCompartments.concat(all_chrom_seg)

	def _do_segment_chromosome(self, chrom, s1, s2):
		_logger = logging.getLogger(self.__class__.__name__)
		_logger.info(f"Processing {chrom}")
		return self.segment_chromosome(s1.get_chromosomes(chrom), s2.get_chromosomes(chrom))


	def segment(self, s1: CalderSubCompartments, s2: CalderSubCompartments,
				n_cpus: Optional[int] = None) -> CalderDifferentialCompartments:
		_logger = logging.getLogger(self.__class__.__name__)
		_logger.info("Binning sample 1")
		s1.binnify(self.binSize)
		_logger.info("Binning sample 2")
		s2.binnify(self.binSize)
		chroms = set(s1.domains.chr).intersection(set(s2.domains.chr))
		n_cpus = n_cpus if n_cpus is not None else self.n_cpus

		with Pool(processes=n_cpus) as pool:
			results = [pool.apply_async(self._do_segment_chromosome, args=(chrom, s1, s2)) for chrom in chroms]
		return CalderDifferentialCompartments.concat([r.get() for r in results])


class CalderRecursiveDifferentialSegmentator(CalderDifferentialSegmentator):
	def __init__(self, binSize: int, statistical_test: str = "gamma", n_cpus: int = 1):
		self._binSize = binSize
		self._stat_test = statistical_test
		self._null_dist = None
		self._n_cpus = n_cpus

	@property
	def binSize(self) -> int:
		return self._binSize

	@property
	def n_cpus(self) -> int:
		return self._n_cpus

	@property
	def control_distribution(self):
		return self._null_dist

	def build_control_distribution(self, pairs: List[Tuple[CalderSubCompartments, CalderSubCompartments]]):

		def __gamma_segmentation(pairs):
			control_dist = []
			for comp1, comp2 in pairs:
				control_dist.append(self.segment(comp1, comp2))
			control_dist = pd.concat([c.segmentation for c in control_dist], axis=0, ignore_index=True)
			control_dist['length'] = control_dist.end - control_dist.start
			control_dist['abs_value'] = control_dist['value'].abs()
			control_dist['n_bins'] = control_dist['length']//self.binSize
			control_dist['statistic'] = control_dist['n_bins']*control_dist['abs_value']
			control_dist = control_dist.dropna(subset=['statistic']).groupby("chr").apply(lambda x: x.statistic.values).to_dict()
			return control_dist

		def __delta_ranks(pairs):
			all_controls = []
			for c1, c2 in pairs:
				c12_rank = self.get_binned_deltaRank(c1, c2)
				all_controls.append(c12_rank)
			all_controls = pd.concat(all_controls, axis=0, ignore_index=True)
			all_controls['abs_delta_rank'] = all_controls['delta_rank'].abs()
			null_dist = all_controls.dropna(subset=['abs_delta_rank']).groupby('chr').apply(lambda x: x.abs_delta_rank.values).to_dict()
			return null_dist

		if self._stat_test == "gamma":
			self._null_dist = __gamma_segmentation(pairs)
		elif self._stat_test == "delta":
			self._null_dist = __delta_ranks(pairs)
		else:
			raise ValueError("Unknown statistical test")
		return self._null_dist

	def get_pvalue(self, core: dict) -> Tuple[float, float]:

		def __gamma_pvalue(chrom, value, n_bins):
			statistic = np.abs(value)*n_bins
			pvalue = ((self._null_dist[chrom] > statistic).sum() + 1)/(self._null_dist[chrom].shape[0] + 1)
			return statistic, pvalue

		def __delta_pvalue(chrom, v):
			if (~np.isnan(v)).sum() > 0:
				statistic = np.max(np.abs( v[~np.isnan(v)] ))
				pvalue = ((self._null_dist[chrom] > statistic).sum() + 1)/(self._null_dist[chrom].shape[0] + 1)
			else:
				statistic = np.NaN
				pvalue = np.NaN
			return statistic, pvalue


		if self._null_dist is None:
			return np.NaN, np.NaN
		else:
			if self._stat_test == "gamma":
				statistic, pvalue = __gamma_pvalue(chrom = core['chromosome'], value = core['value'], n_bins=core['n_bins'])
			elif self._stat_test == "delta":
				statistic, pvalue = __delta_pvalue(chrom = core['chromosome'], v = core['binned_delta_ranks'])
			else:
				raise ValueError("Unknown statistical test")
			return statistic, pvalue


	def segment_chromosome(self,
						   s1: CalderSubCompartments,
						   s2: CalderSubCompartments) -> CalderDifferentialCompartments:


		s12_rank = self.get_binned_deltaRank(s1, s2)
		values = s12_rank["delta_rank"].values

		expected_std = {1: np.inf, values.shape[0]: 0}
		def __get_expected_std(v, n_points=100):
			if v.shape[0] in expected_std.keys():
				return expected_std[v.shape[0]]
			else:
				stds = np.zeros(n_points, dtype=float)
				for p in range(n_points):
					start_p = np.random.choice(values.shape[0] - v.shape[0])
					end_p = start_p + v.shape[0]
					seg = values[start_p : end_p]
					if np.sum(~np.isnan(seg)) > 0:
						std_p = np.nanstd(seg)
						stds[p] = std_p
				expected_std[v.shape[0]] = np.nanmean(stds)
				return expected_std[v.shape[0]]


		def __pos_neg_segmentation(v, segment_value_function = np.mean):
			""" Given a signal, it segments it based on positive and negative values
			It accepts:
			- v: Numpy array
			- segment_value_function: function accepting a numpy array and returning the imputed
				value for each segment, on the basis of the array values
			-
			"""

			segments = []
			seg_values = []
			seg_start = 0
			prev_v = np.NaN
			for i in range(v.shape[0]):
				if not np.isnan(prev_v):
					if np.sign(v[i]) != np.sign(prev_v):
						segments.append([seg_start, i])
						seg_values.append(segment_value_function(v[seg_start:i]))
						seg_start = i
				else:
					seg_start = i
				prev_v = v[i]
			if not np.isnan(prev_v):
				segments.append([seg_start, i + 1])
				seg_values.append(segment_value_function(v[seg_start:i + 1]))
			segments = np.array(segments, dtype=int)
			seg_values = np.array(seg_values)
			return segments, seg_values

		def __recursive_segmentation(X, sub_mean=True):
			X = X.sort_values(['chr', 'start', 'end']).reset_index(drop=True)
			v = X["delta_rank"].values
			rank1 = np.nanmean(X['domain_rank_1'].values)
			rank2 = np.nanmean(X['domain_rank_2'].values)
			mean_signal = np.nanmean(v)
			std_signal = np.nanstd(v)
			exp_std = __get_expected_std(v)
			if std_signal >= exp_std:
				segments, _ = __pos_neg_segmentation(v - mean_signal if sub_mean else v)
				all_rrs = []
				for s in range(segments.shape[0]):
					seg = segments[s, :]
					Xr = X.iloc[seg[0]:seg[1]]
					rrs = __recursive_segmentation(Xr)
					all_rrs.append(rrs)
				result = pd.concat(all_rrs, axis=0, ignore_index=True)
			else:
				result = BedTool.from_dataframe(X.dropna(subset=["delta_rank"])).merge().to_dataframe(names = ['chr', 'start', 'end'])
				result['value'] = mean_signal
				result['rank1'] = rank1
				result['comp1'] = s1.discretize_rank(rank1, chrom = result.chr.iloc[0])
				result['rank2'] = rank2
				result['comp2'] = s2.discretize_rank(rank2, chrom = result.chr.iloc[0])
				statistic, pvalue = self.get_pvalue(core = {
														"chromosome": X.chr.iloc[0],
														"binned_delta_ranks": v,
														"value": mean_signal,
														"length": result.end.max() - result.start.min(),
														"n_bins": (result.end.max() - result.start.min())//self.binSize
													})
				result['statistic'] = statistic
				result['pvalue'] = pvalue
			result = result.sort_values(['chr', 'start', 'end'])
			return result

		result = __recursive_segmentation(s12_rank, sub_mean=False)
		expected_std = [ {"chr": s12_rank.chr.iloc[0], "segment_size": segment_size, "expected_std": exp_std } for segment_size, exp_std in expected_std.items() ]
		expected_std = pd.DataFrame.from_dict(expected_std)
		expected_std = expected_std.sort_values(['chr', "segment_size"]).reset_index(drop=True)
		return CalderDifferentialCompartments(input_df = result,
											  signal_df = s12_rank[["chr", "start", "end",
																   "domain_rank_1", 'domain_rank_2',
																   'compartment_label_8_1', 'compartment_label_8_2',
																   "delta_rank"]],
											  expected_std=expected_std)



segmentator = CalderRecursiveDifferentialSegmentator(binSize = BINSIZE, n_cpus=25)
# segmentator.build_control_distribution([(s1.get_chromosomes('chr1'), s2.get_chromosomes('chr1')) for s1, s2 in control_comparisons])
segmentator.build_control_distribution(control_comparisons)
segments = segmentator.segment(comps[SAMPLE1], comps[SAMPLE2])




segmentator_Delta = CalderRecursiveDifferentialSegmentator(binSize = BINSIZE, statistical_test="delta")
segmentator_Delta.build_control_distribution([(s1.get_chromosomes('chr1'), s2.get_chromosomes('chr1')) for s1, s2 in control_comparisons])
segments_Delta = segmentator_Delta.segment(comps[SAMPLE1].get_chromosomes('chr1'), comps[SAMPLE2].get_chromosomes('chr1'))






sns.scatterplot(data=segments.segmentation.assign(
						n_bins = lambda x: (x.end - x.start)//BINSIZE,
						negLog10Pvalue = lambda x: -1*np.log10(x.pvalue)
					),
				x = "value", y = "negLog10Pvalue", size = "n_bins")
plt.show()



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
shape, loc, scale = gamma.fit(null_dist[null_dist > 0])

# Create a histogram of the data
plt.hist(null_dist[null_dist > 0], bins=50, alpha=0.5, density=True)

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



# COMPARING CBS WITH DIFFCOMP
fig, ax = plt.subplots(1, 1, figsize=(20, 3))

cmap = colormaps['coolwarm']
X = segments.signal.sort_values(['chr', 'start', 'end']).reset_index(drop=True)
v = X["delta_rank"].values
vs  = np.convolve(v, np.ones(4)/4, mode='same')

ax.plot(X.start.values, v, '.', color = '#D78521', markersize=2)
for t in segments.segmentation[segments.segmentation.pvalue <= 0.01].itertuples():
    ax.axvspan(t.start - BINSIZE//2, t.end - BINSIZE//2, alpha=0.2, color = cmap(np.clip((t.value*5 + 1)/2, 0, 1)))
ax.axhline(0, color='black', linewidth=1, linestyle='--')
ax.set_ylabel("$\Delta rank$")
ax.set_title("DiffComp algorithm", loc = "left", size=15)
ax.set_xlabel(f"Genomic coordinate (bp)")
plt.show()
