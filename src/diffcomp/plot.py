import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .diff import CalderDifferentialCompartments

CORE_PALETTE = {"Inactivation": "#0000FF", "Activation": "#FF0000", "NONE": "#808080"}


def __find_pvalue_col(cores: CalderDifferentialCompartments):
    pvalue_col = [c for c in cores.segmentation.columns if c.endswith("_pvalue")]
    if len(pvalue_col) == 0:
        raise ValueError("Missing p-value information")
    else:
        pvalue_col = pvalue_col[0]
    return pvalue_col


def volcanoPlot(
    cores: CalderDifferentialCompartments,
    value_col="value",
    pvalue_col=None,
    size=None,
    absValue_thresh=0,
    pvalue_thresh=0.05,
    pvalue_negLog=True,
    xlim=(-1, 1),
    core_type=None,
    palette=None,
    legend=True,
    title=None,
    ax=None,
    **kwargs,
):
    if pvalue_col is None:
        pvalue_col = __find_pvalue_col(cores)
    pvalue_type = pvalue_col.split("_")[0]

    X = cores.filter(value_col, absValue_thresh, pvalue_col, pvalue_thresh)

    if pvalue_negLog:
        X = X.assign(pvalue=lambda x: -1 * np.log10(x[pvalue_col]))
        pvalue_name = f"$-\log_{{10}}$(p-value) [{pvalue_type}]"
    else:
        X = X.assign(pvalue=x[pvalue_col])
        pvalue_name = f"p-value [{pvalue_type}]"

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        return_fig = True

    palette = CORE_PALETTE if palette is None else palette

    core_type = "core_type" if core_type is None else core_type
    size = "core_size" if size is None else size
    sizes = (5, 50) if size == "core_size" else None
    sns.scatterplot(
        data=X,
        x=value_col,
        y="pvalue",
        hue=core_type,
        size=size,
        sizes=sizes,
        size_norm=(50000, 5000000),
        palette=palette,
        ax=ax,
        legend=legend,
        **kwargs,
    )
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


def transitionHeatmap(
    cores: CalderDifferentialCompartments,
    value_col="value",
    absValue_thresh=0,
    pvalue_col=None,
    pvalue_thresh=0.05,
    ax=None,
    cmap="Greys",
    **kwargs,
):
    if pvalue_col is None:
        pvalue_col = __find_pvalue_col(cores)
    X = cores.filter(value_col, absValue_thresh, pvalue_col, pvalue_thresh)
    X = X[X.significance == "PASS"]
    H = pd.pivot_table(
        X.groupby(["comp1", "comp2"]).size().to_frame("n_cores").reset_index(),
        index="comp2",
        columns="comp1",
        values="n_cores",
        fill_value=0,
    )

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        return_fig = True

    sns.heatmap(
        H, ax=ax, cmap=cmap, cbar=False, annot=True, fmt="d", linewidth=1, **kwargs
    )
    ax.set_xlabel("Starting compartment")
    ax.set_ylabel("Ending compartment")
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

    if return_fig:
        return fig, ax


def barPlot(
    cores: CalderDifferentialCompartments,
    value_col="value",
    absValue_thresh=0,
    pvalue_col=None,
    pvalue_thresh=0.05,
    stat="count_type",
    include_not_significant=False,
    palette=None,
    ax=None,
    **kwargs,
):
    if pvalue_col is None:
        pvalue_col = __find_pvalue_col(cores)

    X = cores.filter(value_col, absValue_thresh, pvalue_col, pvalue_thresh)

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(5, 1.5))
        return_fig = True

    order = ["Activation", "Inactivation"]
    if include_not_significant:
        order.append("NONE")

    palette = CORE_PALETTE if palette is None else palette

    if stat == "count_type":
        sns.countplot(
            data=X, y="core_type", palette=palette, ax=ax, order=order, **kwargs
        )
        ax.set_ylabel("")
        ax.set_xlabel("N. CoREs")
    elif stat == "covered_genome_type":
        sns.barplot(
            data=X,
            x="core_size",
            y="core_type",
            estimator="sum",
            palette=palette,
            ax=ax,
            order=order,
            **kwargs,
        )
        ax.set_ylabel("")
        ax.set_xlabel("Covered genome (bp)")
    ax.set_title("CoREs")
    sns.despine(ax=ax, trim=True)

    if return_fig:
        return fig, ax
