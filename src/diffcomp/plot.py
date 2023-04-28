from .diff import CalderDifferentialCompartments
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


def __to_core_type(x):
	if x.significance == "PASS":
		if x.value < 0:
			return "Inactivation"
		elif x.value > 0:
			return "Activation"
		else:
			raise ValueError("Impossible to have zero value for significant CORE")
	else:
		return "NONE"

CORE_PALETTE = {
	"Inactivation": '#0000FF',
	'Activation': '#FF0000',
	'NONE': '#808080'
}


def volcanoPlot(cores: CalderDifferentialCompartments,
				value_col = "value",
				pvalue_col = None,
				size = None,
				pvalue_thresh=0.05,
				pvalue_negLog = True,
				xlim=(-1, 1),
				core_type = None,
				palette = None,
				legend=True,
				title = None,
				ax = None,
				**kwargs):

	if pvalue_col is None:
		pvalue_col = [c for c in cores.segmentation.columns if c.endswith("_pvalue")]
		if len(pvalue_col) == 0:
			raise ValueError("Missing p-value information")
		else:
			pvalue_col = pvalue_col[0]
	pvalue_type = pvalue_col.split("_")[0]
	if pvalue_negLog:
		X = cores.segmentation.assign(pvalue = lambda x: -1*np.log10(x[pvalue_col]) )
		pvalue_name = f"$-\log_{{10}}$(p-value) [{pvalue_type}]"
	else:
		X = cores.segmentation.assign(pvalue = x[pvalue_col] )
		pvalue_name = f"p-value [{pvalue_type}]"
	if ax is None:
		fig, ax = plt.subplots(1, 1, figsize=(5, 5))
		return_fig = True

	X = X.assign(significance = lambda x: (x[pvalue_col] <= pvalue_thresh).map(lambda y: "PASS" if y else "FAIL"),
				 core_type = lambda x: x.apply(__to_core_type, axis=1),
				 core_size = lambda x: x.end - x.start)

	palette = CORE_PALETTE if palette is None else palette

	core_type = "core_type" if core_type is None else core_type
	size = "core_size" if size is None else size
	sizes = (5, 50) if size == 'core_size' else None
	sns.scatterplot(data = X,
					x = value_col,
					y = "pvalue",
					hue = core_type,
					size = size,
					sizes = sizes,
					size_norm=(50000, 5000000),
					palette=palette,
					ax=ax,
					legend = legend,
					**kwargs)
	ax.set_ylabel(pvalue_name)
	value_name = value_col if value_col != "value" else "Repositioning score"
	ax.set_xlabel(value_name)
	ax.set_xlim(xlim)
	if legend:
		handles, labels = ax.get_legend_handles_labels()
		nh, nl = [], []
		for h, l in zip(handles, labels):
			if l == size:
				break
			nh.append(h)
			nl.append(l)
		ax.legend(nh, nl, bbox_to_anchor=(1, 1))


	if title is not None:
		ax.set_title(title)

	sns.despine(trim=True, ax=ax)

	if return_fig:
		return fig, ax
