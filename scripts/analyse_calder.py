import pandas as pd
from diffcomp.calder import CalderSubCompartments
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

CHROMS = [f'chr{x}' for x in range(1, 23)] + ['chrX']

calder_subcompartments_path = "tests/data/calder/RPE_TP53_Ctrl_mega_CompartmentDomains.bed"

calder_subcompartments = CalderSubCompartments(calder_subcompartments_path)

calder_subcompartments.domains['domain_rank_id'] = calder_subcompartments.domains.groupby("chr").domain_rank.rank().astype(int)


sns.relplot(data = calder_subcompartments.domains,
			col = "chr",
			col_order = CHROMS,
			col_wrap=8,
			x = "domain_rank_id",
			y = "domain_rank",
			facet_kws=dict(sharex=False))
plt.show()

sns.displot(data = calder_subcompartments.domains,
			col = "chr",
			col_order = ['chr1'],
			x = "domain_rank",
			kind = 'hist',
			facet_kws=dict(sharey=False))
plt.show()
