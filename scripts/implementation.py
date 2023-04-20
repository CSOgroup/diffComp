from diffcomp.calder import CalderSubCompartments
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

BINSIZE = 50000


# Compartments
control1 = CalderSubCompartments("tests/data/calder/RPE_TP53_Ctrl_1_CompartmentDomains.bed")
control2 = CalderSubCompartments("tests/data/calder/RPE_TP53_Ctrl_2_CompartmentDomains.bed")

sample1 = CalderSubCompartments("tests/data/calder/RPE_TP53_Ctrl_mega_CompartmentDomains.bed")
sample2 = CalderSubCompartments("tests/data/calder/RPE_TP53_20w0T1_mega_CompartmentDomains.bed")


# Binned compartments
control1_binned = control1.binnify(BINSIZE)
control2_binned = control2.binnify(BINSIZE)
sample1_binned = sample1.binnify(BINSIZE)
sample2_binned = sample2.binnify(BINSIZE)



aligned = sample1_binned.rename(columns = {"compartment_label_8": "comp1", "domain_rank": "rank1"}).drop('compartment_label', axis=1)\
				.merge(sample2_binned.rename(columns = {"compartment_label_8": "comp2", "domain_rank": "rank2"}).drop('compartment_label', axis=1),
						on = ['chr', 'start', 'end'])\
				.assign(deltaRank = lambda x: x.rank2 - x.rank1)

aligned_control = control1_binned.rename(columns = {"compartment_label_8": "comp1", "domain_rank": "rank1"}).drop('compartment_label', axis=1)\
				.merge(control2_binned.rename(columns = {"compartment_label_8": "comp2", "domain_rank": "rank2"}).drop('compartment_label', axis=1),
						on = ['chr', 'start', 'end'])\
				.assign(deltaRank = lambda x: x.rank2 - x.rank1)


sns.displot(data = pd.concat([aligned.assign(name = "observed"), aligned_control.assign(name = "control")], axis=0, ignore_index=True),
			x = "log2RatioRank", hue = "name", kind='kde', common_norm=False)
plt.show()


deltaRanks = aligned.loc[aligned.chr == 'chr1', 'deltaRank'].values
deltaRanks_control = aligned_control.loc[aligned_control.chr == 'chr1', 'deltaRank'].values

changePoints = np.where(np.diff(np.sign(np.nan_to_num(deltaRanks)), prepend=0) != 0)[0]
changePoints_control = np.where(np.diff(np.sign(np.nan_to_num(deltaRanks_control)), prepend=0) != 0)[0]





fig, ax = plt.subplots(2, 1)
ax[0].plot(deltaRanks[:500], marker='.')
ax[0].axhline(0, color='black', linestyle='--', linewidth=1)
[ax[0].axvline(x, color = "red", linewidth=0.1) for x in changePoints[changePoints < 500]]
ax[1].plot(deltaRanks_control[:500], marker='.')
ax[1].axhline(0, color='black', linestyle='--', linewidth=1)
[ax[1].axvline(x, color = "red", linewidth=0.1) for x in changePoints_control[changePoints_control < 500]]
plt.show()
