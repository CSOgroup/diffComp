import pandas as pd
import numpy as np
from typing import Optional
from pybedtools.bedtool import BedTool


LABELS_8_LEVELS = ["{}.{}.{}".format(l1, l2, l3) for l1 in ['A', 'B'] for l2 in [1, 2] for l3 in [1, 2]]
CALDER_SUBCOMPARTMENT_HEADER = ['chr', 'start', 'end', 'compartment_label', 'domain_rank_8_levels', 
								'strand', 'tickStart', 'tickEnd', 'color']
COMPS_HEADER = ['chr', 'start', 'end', 'compartment_label', 'compartment_label_8', 'domain_rank']


class CalderSubCompartments:
	def __init__(self, 
				 path: str = None, 
				 input_df: Optional[pd.DataFrame] = None, 
				 genome: str = "hg19", 
				 binned: Optional[dict] = None,
				 **kwargs):
		if (path is not None) and (input_df is None):
			self._comps = self.read(path, **kwargs)
		elif (path is None) and (input_df is not None):
			if input_df.columns.tolist() == COMPS_HEADER:
				self._comps = input_df
			else:
				raise ValueError("Input dataframe has improper header")
		else:
			raise ValueError("'path' and 'input' are mutually exclusive, and at least one has to be provided")

		self._genome = genome
		if binned is None:
			self._binned = dict()
		else:
			self._binned = binned

	@staticmethod
	def read(path: str, coordinates: str = "zero-based") -> pd.DataFrame:
		""" 
		"""
		s_comps = pd.read_csv(path, sep='\t', header=None,
							  names = CALDER_SUBCOMPARTMENT_HEADER)
		s_comps['compartment_label_8'] = s_comps.compartment_label.map(lambda x: x[:5])

		domain_rank = s_comps.groupby("chr").compartment_label.rank(ascending=False, pct=True)
		s_comps['domain_rank'] = domain_rank

		min_max_domain_rank = s_comps.groupby('chr')['domain_rank']\
									 .min()\
									 .to_frame("min_domain_rank")\
									 .reset_index()\
									 .merge(s_comps.groupby('chr')['domain_rank'].max().to_frame("max_domain_rank").reset_index())

		s_comps = s_comps.merge(min_max_domain_rank)
		s_comps['domain_rank'] = (s_comps.domain_rank - s_comps.min_domain_rank)/(s_comps.max_domain_rank - s_comps.min_domain_rank)

		if coordinates != 'zero-based':
			s_comps['start'] -= 1
			s_comps['tickStart'] -= 1
		s_comps['compartment_label_8'] = pd.Categorical(s_comps['compartment_label_8'], categories=LABELS_8_LEVELS, ordered=True)
		s_comps = s_comps[COMPS_HEADER]
		return s_comps

	@property
	def domains(self) -> pd.DataFrame:
		return self._comps

	def binnify(self, binSize: int):
		if binSize in self._binned.keys():
			return self._binned[binSize]
		else:
			bins = BedTool().window_maker(w=binSize, genome=self._genome)
			comps_binned = bins.map(BedTool.from_dataframe(self._comps).sort(), 
									c=[4, 5, 6], 
									o=['distinct', 'distinct', 'max'])\
				 			   .to_dataframe(names = COMPS_HEADER)
			comps_binned['compartment_label_8'] = comps_binned['compartment_label_8'].map(lambda x: "Unknown" if x == '.' else x)
			comps_binned['compartment_label_8'] = pd.Categorical(comps_binned['compartment_label_8'], categories=LABELS_8_LEVELS, ordered=True)
			comps_binned['compartment_label'] = comps_binned['compartment_label'].map(lambda x: "Unknown" if x == '.' else x)
			comps_binned['domain_rank'] = comps_binned['domain_rank'].map(lambda x: np.NaN if x == '.' else float(x))
			comps_binned = comps_binned[comps_binned.chr.isin(self._comps.chr.unique())].reset_index(drop=True)
			self._binned[binSize] = comps_binned
			return self._binned[binSize]

	def get_chromosome(self, chrom: str):
		return CalderSubCompartments(input_df=self._comps[self._comps.chr == chrom].reset_index(drop=True), 
									 genome = self._genome,
									 binned = {k:v[v.chr == chrom] for k,v in self._binned.items()})




