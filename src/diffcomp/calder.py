from typing import Dict, List, Optional, Union

import bioframe
import numpy as np
import pandas as pd

LABELS_8_LEVELS = [
    "{}.{}.{}".format(l1, l2, l3) for l1 in ["A", "B"] for l2 in [1, 2] for l3 in [1, 2]
]
CALDER_SUBCOMPARTMENT_HEADER = [
    "chr",
    "start",
    "end",
    "compartment_label",
    "domain_rank_8_levels",
    "strand",
    "tickStart",
    "tickEnd",
    "color",
]
COMPS_HEADER = [
    "chr",
    "start",
    "end",
    "compartment_label",
    "compartment_label_8",
    "domain_rank",
]

L8_MIN = [0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875]
L8_MAX = [0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1]
L8_LABS = ["B.2.2", "B.2.1", "B.1.2", "B.1.1", "A.2.2", "A.2.1", "A.1.2", "A.1.1"]

REF_DF = pd.DataFrame(
    {"compartment_label_8": L8_LABS, "new_min": L8_MIN, "new_max": L8_MAX}
)


# rank_correction = True

class CalderSubCompartments:
    def __init__(
        self,
        path: str = None,
        input_df: Optional[pd.DataFrame] = None,
        genome: str = "hg19",
        binned: Optional[dict] = None,
        rank_correction: bool = False,
        **kwargs,
    ):
        if (path is not None) and (input_df is None):
            self._comps = self.read(path, rank_correction, **kwargs)
        elif (path is None) and (input_df is not None):
            if input_df.columns.tolist() == COMPS_HEADER:
                self._comps = input_df
            else:
                raise ValueError("Input dataframe has improper header")
        else:
            raise ValueError(
                "'path' and 'input' are mutually exclusive, and at least one has to be provided"
            )

        self._genome = genome
        if binned is None:
            self._binned = dict()
        else:
            self._binned = binned

        self._rank_ranges = (
            self._comps.groupby(["chr", "compartment_label_8"], observed=False)["domain_rank"] # Set to fix futureWarning
            .min()
            .to_frame("min_rank")
            .reset_index()
        )
        res = []
        for chrom, df in self._rank_ranges.groupby("chr", observed=False):
            df = df.sort_values("compartment_label_8")
            df["max_rank"] = df.min_rank.shift().fillna(1)
            res.append(df)
        res = pd.concat(res, axis=0, ignore_index=True)
        self._rank_ranges = res

    @staticmethod
    def read(path: str, rank_correction: bool = False, coordinates: str = "zero-based") -> pd.DataFrame:
        """ """
        s_comps = pd.read_csv(
            path, sep="\t", header=None, names=CALDER_SUBCOMPARTMENT_HEADER
        )
        s_comps["compartment_label_8"] = s_comps.compartment_label.map(lambda x: x[:5])

        domain_rank = s_comps.groupby("chr").compartment_label.rank(
            ascending=False, pct=True
        )

        s_comps["domain_rank"] = domain_rank

        if rank_correction:
            # Here we normalize to expected intervals to obtain homogeneous ranks
            # across different trees

            # Add a message to stdOut
            print(
                "Rank correction is enabled. Normalizing ranks to obtain homogeneous intervals"
            )

            min_max = (
                s_comps.groupby(["chr", "compartment_label_8"])
                .agg(old_min=("domain_rank", "min"), old_max=("domain_rank", "max"))
                .reset_index()
            )

            # s_comps = s_comps.merge(min_max_domain_rank)
            s_comps = s_comps.merge(min_max)
            s_comps = s_comps.merge(REF_DF)

            # new_CompRank = new_min + ((comp_rank-old_min)/(old_max-old_min))*(new_max-new_min)

            s_comps["domain_rank"] = s_comps.new_min + (
                (s_comps.domain_rank - s_comps.old_min)
                / (s_comps.old_max - s_comps.old_min)
            ) * (s_comps.new_max - s_comps.new_min)

        else:
            min_max_domain_rank = (
                s_comps.groupby("chr")["domain_rank"]
                .min()
                .to_frame("min_domain_rank")
                .reset_index()
                .merge(
                    s_comps.groupby("chr")["domain_rank"]
                    .max()
                    .to_frame("max_domain_rank")
                    .reset_index()
                )
            )

            s_comps = s_comps.merge(min_max_domain_rank)
            s_comps["domain_rank"] = (s_comps.domain_rank - s_comps.min_domain_rank) / (
                s_comps.max_domain_rank - s_comps.min_domain_rank
            )

        if coordinates != "zero-based":
            s_comps["start"] -= 1
            s_comps["tickStart"] -= 1
        s_comps["compartment_label_8"] = pd.Categorical(
            s_comps["compartment_label_8"], categories=LABELS_8_LEVELS, ordered=True
        )
        s_comps = s_comps[COMPS_HEADER]
        return s_comps

    @property
    def chromosomes(self) -> List[str]:
        return self._comps.chr.unique().tolist()

    @property
    def domains(self) -> pd.DataFrame:
        return self._comps

    def get_domain_boundaries(self, chrom=None) -> Union[Dict[str, list], list]:
        if chrom is None:
            res = {}
            for c, df in self._comps.groupby("chr"):
                c_bounds = np.array(
                    sorted(set(df["start"].tolist() + df["end"].tolist()))
                )
                res[c] = c_bounds
        else:
            df = self._comps[self._comps.chr == chrom]
            res = np.array(sorted(set(df["start"].tolist() + df["end"].tolist())))
        return res

    def binnify(self, binSize: int):
        if binSize in self._binned.keys():
            return self._binned[binSize]
        else:
            bins = bioframe.binnify(
                bioframe.fetch_chromsizes(self._genome), binsize=binSize
            )
            comps_binned = (
                bioframe.overlap(bins, self._comps.rename(columns={"chr": "chrom"}))
                .drop(["chrom_", "start_", "end_"], axis=1)
                .rename(columns=lambda x: x[:-1] if x.endswith("_") else x)
                .rename(columns={"chrom": "chr"})
                .assign(
                    compartment_label=lambda y: y.compartment_label.fillna("Unknown")
                )
            )
            comps_binned = comps_binned[
                comps_binned.chr.isin(self._comps.chr.unique())
            ].reset_index(drop=True)
            self._binned[binSize] = comps_binned
            return self._binned[binSize]

    def get_chromosomes(self, chrom: Union[str, List[str]]):
        if isinstance(chrom, str):
            return CalderSubCompartments(
                input_df=self._comps[self._comps.chr == chrom].reset_index(drop=True),
                genome=self._genome,
                binned={k: v[v.chr == chrom] for k, v in self._binned.items()},
            )
        elif isinstance(chrom, list):
            return CalderSubCompartments(
                input_df=self._comps[self._comps.chr.isin(chrom)].reset_index(
                    drop=True
                ),
                genome=self._genome,
                binned={k: v[v.chr.isin(chrom)] for k, v in self._binned.items()},
            )
        else:
            raise ValueError("Unknown chrom type")

    def discretize_rank(self, rank: float, chrom: Optional[str] = None):
        if chrom is None:
            intervals = self._rank_ranges[
                (self._rank_ranges.min_rank <= rank)
                & (rank <= self._rank_ranges.max_rank)
            ]
        else:
            intervals = self._rank_ranges[
                (self._rank_ranges.chr == chrom)
                & (self._rank_ranges.min_rank <= rank)
                & (rank <= self._rank_ranges.max_rank)
            ]
        return intervals.compartment_label_8.mode().iloc[0]
