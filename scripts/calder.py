from cooler import Cooler
import numpy as np

import matplotlib.pyplot as plt


def load_chrom_ints_as_sparse_matrix(cool, chrom, balance=None):
	M = cool.matrix(balance=balance).fetch(chrom)
	return M


def columnwise_pearson_correlation(M_log):
	# Scaling the columns
	# scaler = StandardScaler()
	# M_scaled = scaler.fit_transform(M_log)
	M_scaled = np.nan_to_num((M_log - np.nanmean(M_log, axis=0))/np.nanstd(M_log, axis=0))
	# M_scaled = np.ma.array(M_scaled, mask = np.isnan(M_scaled))
	P = 1/(M_scaled.shape[0])* M_scaled.T.dot(M_scaled)
	return P

binSize=50000
balance = "KR"
chrom = "10"
cool_path=f"data/calder/inter_30.mcool::/resolutions/{binSize}"

cool = Cooler(cool_path)
M = load_chrom_ints_as_sparse_matrix(cool, chrom, balance=balance)

# Taking the log of the data
M_log = np.log2(M + 1)

plt.matshow(M_log, cmap="Reds")
plt.title("Log2(contacts)")


# Computing pearson correlation for each column
M_corr = columnwise_pearson_correlation(M_log)

plt.matshow(M_corr, cmap='bwr', vmin=-1, vmax=1)
plt.title("Corr(Log2(contacts))")

M_corr_corr = columnwise_pearson_correlation(M_corr)

plt.matshow(M_corr_corr, cmap='bwr', vmin=-1, vmax=1)
plt.title("Corr(Corr(Log2(contacts)))")


M_atah = np.arctanh(M_corr_corr)

plt.matshow(M_atah, cmap='bwr', vmin=-2, vmax=2)
plt.title("arctanh(Corr(Corr(Log2(contacts))))")
plt.show()
