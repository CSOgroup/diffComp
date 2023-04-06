# -*- coding: utf-8 -*-
# @Author: lucananni93
# @Date:   2022-04-29 13:24:38
# @Last Modified by:   lucananni93
# @Last Modified time: 2022-04-29 18:40:25

import ruptures as rpt
from diffcomp.calder import CalderSubCompartments
from diffcomp.diff import CalderRecursiveDifferentialSegmentator
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

FILE1="tests/data/calder/RPE_TP53_Ctrl_mega_CompartmentDomains.bed"
FILE2="tests/data/calder/RPE_TP53_20w0T1_mega_CompartmentDomains.bed"
CHROM = "chr5"

comps1 = CalderSubCompartments(FILE1).get_chromosome(CHROM)
comps2 = CalderSubCompartments(FILE2).get_chromosome(CHROM)

segmentator = CalderRecursiveDifferentialSegmentator(50000, 0.1)
CoRE_segments = segmentator.segment(comps1, comps2).segmentation

comps12 = segmentator.get_binned_deltaRank(comps1, comps2).dropna(subset=['delta_rank'])


model = rpt.Binseg(model="l2")


def CalderChangePointDifferentialSegmentator(X, model, predict_kws=None):

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


	comps12 = X.copy()

	Algo_segments = []
	comps12['roi'] = (comps12.start != comps12.end.shift()).cumsum()

	for roi, df in comps12.groupby("roi"):
		roi_segments = get_change_points(
								df,
								model,
								predict_kws=predict_kws)
		Algo_segments.append(roi_segments)
	Algo_segments = pd.concat(Algo_segments, axis=0, ignore_index=True)
	return Algo_segments


best_segments = None
best_pen = None
best_distance = None
for pen in np.linspace(0, 0.2, 20):
	Algo_segments = CalderChangePointDifferentialSegmentator(
							comps12,
							model,
							predict_kws=dict(pen = pen))
	print(f"Pen:{pen:<40}N. segments: {Algo_segments.shape[0]}")
	distance = abs(Algo_segments.shape[0] - CoRE_segments.shape[0])
	if (best_segments is None) or (best_distance > distance):
		best_segments = Algo_segments
		best_pen = pen
		best_distance = distance



plt.plot(comps12.start, comps12.delta_rank, 'bo', markersize=2)
[plt.axvline(x - 25000, color = "black", linewidth=1, linestyle = "--") for x in best_segments.start]
[plt.axvline(x - 25000, color = "red", linewidth=1, linestyle = "--") for x in CoRE_segments.start]
plt.show()
