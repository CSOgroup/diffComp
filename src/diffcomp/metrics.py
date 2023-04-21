import numpy as np
import pandas as pd
from abc import ABC, abstractmethod

class SegmentationPairwiseMetric(ABC):
	@staticmethod
	def _check_header(s: pd.DataFrame):
		if s.columns[:3].tolist() != ["chr", 'start', 'end']:
			raise ValueError("Segmentation comparison is available only between dataframes containing 'chr', 'start' and 'end' columns")

	@abstractmethod
	def compare_segments(self, s1: pd.DataFrame, s2: pd.DataFrame):
		pass

	def compare(self,
				s1: pd.DataFrame,
				s2: pd.DataFrame,
				by_chrom=False,
				estimator=None):
		self._check_header(s1)
		self._check_header(s2)
		result = {}
		for chrom in set(s1.chr).intersection(set(s2.chr)):
			chrom_metric = self.compare_segments(s1[s1.chr == chrom], s2[s2.chr == chrom])
			result[chrom] = chrom_metric
		if by_chrom:
			return result
		else:
			vals = list(result.values())
			if estimator is None:
				return np.mean(vals)
			else:
				return estimator(vals)


def moc(cds_p, cds_q):
	"""Measure of concordance between two TAD sets, as defined in https://www.nature.com/articles/s41588-018-0338-y
	If P and Q are two sets of TADs, Np and Nq are the number of TADs in the two sets respectively:

				|-- 1 	if Np = Nq = 1
	MoC(P, Q) = |
				|-- 1 / (SQRT(NpNq) - 1) * (SUM_{i: 1->Np}SUM_{j: 1-> Nq} ( ||Fij||**2 / (||Pi||*||Qj||) )) otherwise

	where Pi and Qj are two individual TADs in P and Q having size ||Pi|| and ||Qj|| respectively (measured in base pairs).
	Finally, ||Fij|| is the size (n. base pairs) of overlap between Pi and Pj.

	This function takes as input two sets of TADs, defined as a list of tuples [(start, end), ...].
	TADs have to belong to the same chromosome, for obvious reasons.
	"""

	def __fill_gaps(tads):
		fills = []
		sorted_tads = sorted(tads, key = lambda x: x[0])
		prev = None
		for t in sorted_tads:
			if prev is not None:
				if t[0] > prev:
					fills.append( (prev, t[0]) )
				elif t[0] < prev:
					raise ValueError(f"Overlapping regions")
			prev = t[1]
		return fills

	def __make_same_limits(tads1, tads2):
		sorted_tads1 = sorted(tads1, key = lambda x: x[0])
		sorted_tads2 = sorted(tads2, key = lambda x: x[0])
		min1, min2 = sorted_tads1[0][0], sorted_tads2[0][0]
		if min1 != min2:
			if min1 < min2:
				sorted_tads2 = [(min1, min2)] + sorted_tads2
			else:
				sorted_tads1 = [(min2, min1)] + sorted_tads1
		max1, max2 = sorted_tads1[-1][1], sorted_tads2[-1][1]
		if max1 != max2:
			if max1 > max2:
				sorted_tads2.append( (max2, max1) )
			else:
				sorted_tads1.append( (max1, max2) )
		return sorted_tads1, sorted_tads2

	def __get_overlap(cd_i, cd_j):
		return max(0, min(cd_i[1], cd_j[1]) - max(cd_i[0], cd_j[0]))

	def __get_size(cd):
		return cd[1] - cd[0]

	def __aggregate_overlaps(cds_p, cds_q):
		r = 0
		for cd_i in cds_p:
			for cd_j in cds_q:
				r += (__get_overlap(cd_i, cd_j)**2)/(__get_size(cd_i)*__get_size(cd_j))
		return r

	cds_p, cds_q = __make_same_limits(cds_p + __fill_gaps(cds_p), cds_q + __fill_gaps(cds_q))
	moc = None
	Np = len(cds_p)
	Nq = len(cds_q)
	if Np == Nq == 1:
		moc = 1
	else:
		moc = (1/(np.sqrt(Np*Nq) - 1))*(__aggregate_overlaps(cds_p, cds_q) - 1)
	return moc


class MeasureOfConcordance(SegmentationPairwiseMetric):

	def compare_segments(self, s1: pd.DataFrame, s2: pd.DataFrame):
		tads1 = list(s1[['start', 'end']].to_records(index=False))
		tads2 = list(s2[['start', 'end']].to_records(index=False))
		return moc(tads1, tads2)
