from __future__ import annotations
from abc import ABC, abstractmethod
from .calder import CalderSubCompartments
from typing import Generic, Optional, List, Tuple
import pandas as pd
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from pybedtools.bedtool import BedTool
import logging
import ruptures as rpt


class CalderDifferentialCompartments:
	CALDER_DIFFERENTIAL_COMPARTMENTS_HEADER = ['chr', 'start', 'end', 'value', "rank1", "comp1", "rank2", "comp2", 'pvalue']

	def __init__(self,
				 path: str = None,
				 input_df: Optional[pd.DataFrame] = None,
				 signal_path: Optional[str] = None,
				 signal_df: Optional[pd.DataFrame] = None):
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
		return CalderDifferentialCompartments(input_df=all_segs, signal_df = all_signals)

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

	def get_binned_deltaRank(self, s1: CalderSubCompartments, s2: CalderSubCompartments) -> pd.DataFrame:
		s1_binned, s2_binned = s1.binnify(self._binSize), s2.binnify(self._binSize)
		s1_rank = s1_binned[['chr', 'start', 'end', 'domain_rank', 'compartment_label_8']]
		s2_rank = s2_binned[['chr', 'start', 'end', 'domain_rank', 'compartment_label_8']]
		s12_rank = s1_rank.merge(s2_rank, on=['chr', 'start', 'end'], suffixes=("_1", "_2"))
		s12_rank['delta_rank'] = s12_rank.domain_rank_2 - s12_rank.domain_rank_1
		return s12_rank

	def segment(self,
				s1: CalderSubCompartments,
				s2: CalderSubCompartments) -> CalderDifferentialCompartments:
		_logger = logging.getLogger(self.__class__.__name__)
		_logger.info("Binning sample 1")
		s1.binnify(self.binSize)
		_logger.info("Binning sample 2")
		s2.binnify(self.binSize)
		chroms = set(s1.domains.chr).intersection(set(s2.domains.chr))
		all_chrom_seg = []
		for chrom in chroms:
			_logger.info(f"Processing {chrom}")
			chrom_seg = self.segment_chromosome(s1.get_chromosomes(chrom), s2.get_chromosomes(chrom))
			all_chrom_seg.append(chrom_seg)
		return CalderDifferentialCompartments.concat(all_chrom_seg)


class CalderChangePointDifferentialSegmentator(CalderDifferentialSegmentator):
	def __init__(self, binSize: int,
				 model: str,
				 model_kws: Optional[dict],
				 predict_kws: Optional[dict]):
		self._binSize = binSize

		if model_kws is None:
			model_kws = dict()
		self._model_kws = model_kws
		if predict_kws is None:
			predict_kws = dict()
		self._predict_kws = predict_kws

		if model == "binary":
			self._model = rpt.Binseg(**self.model_kws)
		else:
			raise ValueError(f"Unknown Change Point detection method ({model})")

	@property
	def binSize(self) -> int:
		return self._binSize

	@staticmethod
	def get_change_points(X, model, predict_kws=None):
		comps12 = X.sort_values(['chr', 'start', 'end']).reset_index(drop=True).copy()
		deltaRank = comps12["delta_rank"].values
		algo = model.fit(deltaRank)
		bkp = algo.predict(**predict_kws)


		comps12['binPos'] = "segment"
		comps12.loc[bkp[:-1], 'binPos'] = 'breakpoint'

		comps12['segment_id'] = (comps12.binPos != comps12.binPos.shift()).cumsum()

		segments = comps12.groupby("segment_id").agg({
									"chr": "first",
									"start": "min",
									"end": "max",
									"delta_rank": "mean"
								})\
								.assign(pvalue = np.NaN)\
								.sort_values(['chr', 'start', 'end'])\
								.reset_index(drop=True)
		return segments

	def segment_chromosome(self,
						   s1: CalderSubCompartments,
						   s2: CalderSubCompartments) -> CalderDifferentialCompartments:
		s12_rank = self.get_binned_deltaRank(s1, s2)
		Algo_segments = []
		s12_rank['roi'] = (s12_rank.start != s12_rank.end.shift()).cumsum()
		for roi, df in s12_rank.groupby("roi"):
			roi_segments = get_change_points(
									df,
									model,
									predict_kws=self.predict_kws)
			Algo_segments.append(roi_segments)
		Algo_segments = pd.concat(Algo_segments, axis=0, ignore_index=True)
		return CalderDifferentialCompartments(input_df = Algo_segments, signal_df = s12_rank[["chr", "start", "end", "delta_rank"]])


class CalderRecursiveDifferentialSegmentator(CalderDifferentialSegmentator):
	def __init__(self, binSize: int, min_std: float = 0.1):
		self._min_std = min_std
		self._binSize = binSize
		self._null_dist = None

	@property
	def binSize(self) -> int:
		return self._binSize

	@property
	def minStd(self) -> float:
		return self._min_std

	@minStd.setter
	def minStd(self, value) -> None:
		self._min_std = value


	def build_control_distribution(self, pairs: List[Tuple[CalderSubCompartments, CalderSubCompartments]]):
		all_controls = []
		for c1, c2 in pairs:
			c12_rank = self.get_binned_deltaRank(c1, c2)
			all_controls.append(c12_rank)
		all_controls = pd.concat(all_controls, axis=0, ignore_index=True)
		null_dist = all_controls.dropna(subset=['delta_rank']).groupby('chr').apply(lambda x: x.delta_rank.values).to_dict()
		self._null_dist = null_dist
		return null_dist


	def segment_chromosome(self,
						   s1: CalderSubCompartments,
						   s2: CalderSubCompartments) -> CalderDifferentialCompartments:
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
			if std_signal > self._min_std:
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

				if (control_dist_abs is not None) and ((~np.isnan(v)).sum() > 0):
					pvalue = ((control_dist_abs[X.chr.iloc[0]] > np.max(np.abs( v[~np.isnan(v)] ))).sum() + 1)/(control_dist_abs[X.chr.iloc[0]].shape[0] + 1)
					result['pvalue'] = pvalue
				else:
					result['pvalue'] = np.NaN
			result = result.sort_values(['chr', 'start', 'end'])
			return result

		s12_rank = self.get_binned_deltaRank(s1, s2)
		control_dist_abs = {chrom:np.abs(v) for chrom, v in self._null_dist.items()} if self._null_dist is not None else None
		result = __recursive_segmentation(s12_rank)
		return CalderDifferentialCompartments(input_df = result, signal_df = s12_rank[["chr", "start", "end",
																					   "domain_rank_1", 'domain_rank_2',
																					   'compartment_label_8_1', 'compartment_label_8_2',
																					   "delta_rank"]])
