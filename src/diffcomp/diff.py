from __future__ import annotations

import logging
import pickle
from abc import ABC, abstractmethod

# import ruptures as rpt
from multiprocessing import Pool
from typing import Generic, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import gamma

from .calder import LABELS_8_LEVELS, CalderSubCompartments


class CalderDifferentialCompartments:
    CALDER_DIFFERENTIAL_COMPARTMENTS_HEADER = [
        "chr",
        "start",
        "end",
        "value",
        "rank1",
        "comp1",
        "rank2",
        "comp2",
    ]

    def __init__(
        self,
        path: str = None,
        input_df: Optional[pd.DataFrame] = None,
        signal_path: Optional[str] = None,
        signal_df: Optional[pd.DataFrame] = None,
        expected_std: Optional[pd.DataFrame] = None,
        chroms: Optional[List[str]] = None,
    ):
        if (path is None) and (input_df is not None):
            self._check_header(input_df)
            self._segs = input_df
        elif (path is not None) and (input_df is None):
            self._segs = self.read(path)
        else:
            raise ValueError(
                "'path' and 'input_df' are mutually exclusive, and at least one has to be provided"
            )

        if chroms is None:
            chroms = self._segs.chr.unique().tolist()
        else:
            if len(set(chroms).difference(set(self._segs.chr.unique().tolist()))) != 0:
                raise ValueError("Specified chromosomes are not present in the dataset")
        self._segs["chr"] = pd.Categorical(
            self._segs["chr"], categories=chroms, ordered=True
        )

        if (signal_path is None) and (signal_df is not None):
            self._signal = signal_df
        elif (signal_path is not None) and (signal_df is None):
            self._signal = self.read_signal(signal_path)
        elif (signal_path is None) and (signal_df is None):
            self._signal = None
        else:
            raise ValueError("'signal_path' and 'signal_df' are mutually exclusive")

        if expected_std is not None:
            self._expected_std = expected_std
        else:
            self._expected_std = None

    @classmethod
    def read(cls, path: str):
        segs = pd.read_csv(path, sep="\t")
        cls._check_header(segs)
        return segs

    @classmethod
    def read_signal(cls, path: str):
        signal = pd.read_csv(path, sep="\t")
        return signal

    @classmethod
    def _check_header(cls, segs: pd.DataFrame):
        if (
            segs.columns[: len(cls.CALDER_DIFFERENTIAL_COMPARTMENTS_HEADER)].tolist()
            != cls.CALDER_DIFFERENTIAL_COMPARTMENTS_HEADER
        ):
            raise ValueError(
                "Calder differential segmentation file must have as initial columns "
                + ", ".join(cls.CALDER_DIFFERENTIAL_COMPARTMENTS_HEADER)
            )
        segs["comp1"] = pd.Categorical(
            segs["comp1"].astype(str), categories=LABELS_8_LEVELS, ordered=True
        )
        segs["comp2"] = pd.Categorical(
            segs["comp2"].astype(str), categories=LABELS_8_LEVELS, ordered=True
        )
        if "core_size" not in segs.columns:
            segs["core_size"] = segs.end - segs.start

    @property
    def segmentation(self):
        return self._segs

    @property
    def signal(self):
        return self._signal

    @property
    def expected_std_by_segment_size(self):
        return self._expected_std

    @staticmethod
    def concat(
        segmentations: List[CalderDifferentialCompartments], **kwargs
    ) -> CalderDifferentialCompartments:
        all_segs = [x.segmentation for x in segmentations]
        all_segs = pd.concat(all_segs, axis=0, ignore_index=True)
        all_segs = all_segs.sort_values(["chr", "start", "end"]).reset_index(drop=True)
        all_signals = [x.signal for x in segmentations if x.signal is not None]
        if len(all_signals) > 0:
            all_signals = pd.concat(all_signals, axis=0, ignore_index=True)
            all_signals = all_signals.sort_values(["chr", "start", "end"]).reset_index(
                drop=True
            )
        else:
            all_signals = None
        all_expected_std = [
            x.expected_std_by_segment_size
            for x in segmentations
            if x.expected_std_by_segment_size is not None
        ]
        if len(all_expected_std) > 0:
            all_expected_std = pd.concat(all_expected_std, axis=0, ignore_index=True)
        else:
            all_expected_std = None
        return CalderDifferentialCompartments(
            input_df=all_segs,
            signal_df=all_signals,
            expected_std=all_expected_std,
            **kwargs,
        )

    def to_tsv(self, path: str):
        self._segs.to_csv(path, sep="\t", index=False, header=True)

    def to_bed(self, path: str):
        cmap = LinearSegmentedColormap.from_list(
            "CoREs_colormap", ["blue", "lightgrey", "red"]
        )
        self._segs.assign(
            name=lambda x: x.chr.astype(str)
            + ":"
            + x.start.astype(str)
            + "-"
            + x.end.astype(str),
            norm_value=lambda x: (x.value + 1) / 2,
            score=lambda x: x.value,
            strand=".",
            thickStart=lambda x: x.start,
            thickEnd=lambda x: x.end,
            color=lambda x: list(
                map(
                    lambda x: ",".join(x),
                    (cmap(x.norm_value.values) * 255)
                    .astype(int)
                    .astype(str)[:, :-1]
                    .tolist(),
                )
            ),
        )[
            [
                "chr",
                "start",
                "end",
                "name",
                "score",
                "strand",
                "thickStart",
                "thickEnd",
                "color",
            ]
        ].to_csv(path, sep="\t", index=False, header=False)

    def write_signal(self, path: str):
        if self._signal is not None:
            self._signal.dropna().to_csv(path, sep="\t", index=False, header=False)

    def get_chromosomes(self, chroms: Union[str, List[chr]]):
        if isinstance(chroms, str):
            chroms = [chroms]
        segs = self._segs[self._segs.chr.isin(chroms)].reset_index(drop=True)
        signal = (
            self._signal[self._signal.chr.isin(chroms)].reset_index(drop=True)
            if self._signal is not None
            else None
        )
        expected_std = (
            self._expected_std[self._expected_std.chr.isin(chroms)].reset_index(
                drop=True
            )
            if self._expected_std is not None
            else None
        )
        return CalderDifferentialCompartments(
            input_df=segs, signal_df=signal, expected_std=expected_std
        )

    def filter(
        self,
        value_col: str = "value",
        absValue_thresh: float = 0,
        pvalue_col: str = "gamma_pvalue",
        pvalue_thresh: float = 0.05,
    ):
        def __to_core_type(x):
            if x.significance == "PASS":
                if x.value < 0:
                    return "Inactivation"
                elif x.value > 0:
                    return "Activation"
                else:
                    raise ValueError(
                        "Impossible to have zero value for significant CORE"
                    )
            else:
                return "NONE"

        X = self._segs.copy()
        X["significance"] = (
            (X[pvalue_col] <= pvalue_thresh) & (X[value_col].abs() >= absValue_thresh)
        ).map(lambda y: "PASS" if y else "FAIL")
        X["core_type"] = X.apply(__to_core_type, axis=1)
        return X

    def __repr__(self):
        result = "------------------------------------\n"
        result += "| Compartment Repositioning Events |\n"
        result += "------------------------------------\n"
        result += f"- N. segments: {self.segmentation.shape[0]}\n\n"
        result += "First segments: \n"
        segmentation_head_str = str(self.segmentation.head())
        segmentation_head_linewidth = len(segmentation_head_str.split("\n")[0])

        result += chr(173) * segmentation_head_linewidth + "\n"
        result += segmentation_head_str + "\n"
        result += ".......\n"
        result += chr(173) * segmentation_head_linewidth + "\n"
        return result


class CalderDifferentialSegmentator(ABC):
    @abstractmethod
    def segment_chromosome(
        self, s1: CalderSubCompartments, s2: CalderSubCompartments
    ) -> CalderDifferentialCompartments:
        pass

    @property
    @abstractmethod
    def binSize(self) -> int:
        pass

    @property
    @abstractmethod
    def n_cpus(self) -> int:
        pass

    def get_binned_deltaRank(
        self, s1: CalderSubCompartments, s2: CalderSubCompartments
    ) -> pd.DataFrame:
        s1_binned, s2_binned = s1.binnify(self._binSize), s2.binnify(self._binSize)
        s1_rank = s1_binned[
            ["chr", "start", "end", "domain_rank", "compartment_label_8"]
        ]
        s2_rank = s2_binned[
            ["chr", "start", "end", "domain_rank", "compartment_label_8"]
        ]
        s12_rank = s1_rank.merge(
            s2_rank, on=["chr", "start", "end"], suffixes=("_1", "_2")
        )
        s12_rank["delta_rank"] = s12_rank.domain_rank_2 - s12_rank.domain_rank_1
        return s12_rank

    # def segment(self,
    #           s1: CalderSubCompartments,
    #           s2: CalderSubCompartments) -> CalderDifferentialCompartments:
    #   _logger = logging.getLogger(self.__class__.__name__)
    #   _logger.info("Binning sample 1")
    #   s1.binnify(self.binSize)
    #   _logger.info("Binning sample 2")
    #   s2.binnify(self.binSize)
    #   chroms = set(s1.domains.chr).intersection(set(s2.domains.chr))
    #   all_chrom_seg = []
    #   for chrom in chroms:
    #       _logger.info(f"Processing {chrom}")
    #       chrom_seg = self.segment_chromosome(s1.get_chromosomes(chrom), s2.get_chromosomes(chrom))
    #       all_chrom_seg.append(chrom_seg)
    #   return CalderDifferentialCompartments.concat(all_chrom_seg)

    def _do_segment_chromosome(self, chrom, s1, s2):
        _logger = logging.getLogger(self.__class__.__name__)
        _logger.info(f"Processing {chrom}")
        return self.segment_chromosome(
            s1.get_chromosomes(chrom), s2.get_chromosomes(chrom)
        )

    def segment(
        self,
        s1: CalderSubCompartments,
        s2: CalderSubCompartments,
        n_cpus: Optional[int] = None,
        chroms: Optional[List[str]] = None,
    ) -> CalderDifferentialCompartments:
        _logger = logging.getLogger(self.__class__.__name__)
        _logger.debug("Binning sample 1")
        s1.binnify(self.binSize)
        _logger.debug("Binning sample 2")
        s2.binnify(self.binSize)
        shared_chroms = set(s1.domains.chr).intersection(set(s2.domains.chr))
        if chroms is None:
            chroms = sorted(shared_chroms)
        else:
            if len(set(chroms).difference(shared_chroms)) != 0:
                raise ValueError(
                    f"Specified chromosomes are not present in the datasets"
                )
        n_cpus = n_cpus if n_cpus is not None else self.n_cpus
        _logger.debug(f"Running segmentation with {n_cpus} CPUs")
        results = [self._do_segment_chromosome(chrom, s1, s2) for chrom in chroms]
        return CalderDifferentialCompartments.concat(results, chroms=chroms)


class CalderRecursiveDifferentialSegmentator(CalderDifferentialSegmentator):
    TESTS = ["gamma", "gamma-fit", "delta"]

    def __init__(
        self,
        binSize: int,
        statistical_test: Union[str, List[str]] = "gamma",
        min_std: Optional[float] = None,
        n_cpus: int = 1,
        random_state: Optional[int] = None,
    ):
        self._binSize = binSize

        if isinstance(statistical_test, str):
            if statistical_test == "all":
                statistical_test = self.TESTS
            elif statistical_test not in self.TESTS:
                raise ValueError("Unknown statistical test")
            else:
                statistical_test = [statistical_test]
        else:
            if len(set(statistical_test).difference(set(self.TESTS))) != 0:
                raise ValueError("Unknown statistical test")

        self._stat_test = statistical_test
        self._null_dist = None
        self._n_cpus = n_cpus
        self._min_std = min_std
        self._logger = logging.getLogger(self.__class__.__name__)
        self._random_generator = np.random.default_rng(random_state)

    @property
    def binSize(self) -> int:
        return self._binSize

    @property
    def n_cpus(self) -> int:
        return self._n_cpus

    @property
    def control_distribution(self):
        return self._null_dist

    def build_control_distribution(
        self, pairs: List[Tuple[CalderSubCompartments, CalderSubCompartments]]
    ):
        res = {}

        def __gamma_empirical(pairs):
            control_dist = []
            for i, pair in enumerate(pairs):
                self._logger.info(f"Control comparison {i + 1}")
                comp1, comp2 = pair
                control_dist.append(self.segment(comp1, comp2))
            control_dist = pd.concat(
                [c.segmentation for c in control_dist], axis=0, ignore_index=True
            )
            control_dist["length"] = control_dist.end - control_dist.start
            control_dist["abs_value"] = control_dist["value"].abs()
            control_dist["n_bins"] = control_dist["length"] // self.binSize
            control_dist["statistic"] = (
                control_dist["n_bins"] * control_dist["abs_value"]
            )
            segments = control_dist.copy()
            control_dist = (
                control_dist.dropna(subset=["statistic"])
                .groupby("chr")
                .apply(lambda x: x.statistic.values)
                .to_dict()
            )
            return control_dist, segments

        def __gamma_fit(pairs):
            if "gamma" in res.keys():
                control_dist = res["gamma"]
            else:
                control_dist, segments = __gamma_empirical(pairs)
                res["gamma"] = control_dist
                res["gamma_segments"] = segments
            control_dist = {
                chrom: gamma.fit(vs[vs > 0]) for chrom, vs in control_dist.items()
            }
            return control_dist

        def __delta_ranks(pairs):
            all_controls = []
            for c1, c2 in pairs:
                c12_rank = self.get_binned_deltaRank(c1, c2)
                all_controls.append(c12_rank)
            all_controls = pd.concat(all_controls, axis=0, ignore_index=True)
            all_controls["abs_delta_rank"] = all_controls["delta_rank"].abs()
            null_dist = (
                all_controls.dropna(subset=["abs_delta_rank"])
                .groupby("chr", observed=False)
                .apply(lambda x: x.abs_delta_rank.values)
                .to_dict()
            )
            return null_dist

        CONTROL_DISTRIBUTIONS = {
            "gamma": __gamma_empirical,
            "gamma-fit": __gamma_fit,
            "delta": __delta_ranks,
        }

        for test in self._stat_test:
            self._logger.info(f"Computing control {test} distribution")
            if test not in res.keys():
                if test == "gamma":
                    cdist, segs = CONTROL_DISTRIBUTIONS[test](pairs)
                    res["gamma"] = cdist
                    res["gamma_segments"] = segs
                else:
                    res[test] = CONTROL_DISTRIBUTIONS[test](pairs)
        self._null_dist = res
        return self._null_dist

    def write_control_distribution(self, path: str):
        if self._null_dist is not None:
            with open(path, "wb") as f:
                pickle.dump(self._null_dist, f)
        else:
            raise ValueError("Control distribution not calculated yet.")

    def load_control_distribution(self, path: str):
        with open(path, "rb") as f:
            null_dist = pickle.load(f)
        # checking if the object is indeed a control distribution
        if not isinstance(null_dist, dict):
            raise ValueError("Control distribution must be a dictionary")
        for test in null_dist.keys():
            if test not in self.TESTS + ["gamma_segments"]:
                raise ValueError(f"{test} is not a valid statistical test")
        self._null_dist = null_dist

    def get_pvalue(self, core: dict) -> Tuple[float, float]:
        def __gamma_pvalue(core):
            chrom, value, n_bins = core["chromosome"], core["value"], core["n_bins"]
            statistic = np.abs(value) * n_bins
            pvalue = ((self._null_dist["gamma"][chrom] > statistic).sum() + 1) / (
                self._null_dist["gamma"][chrom].shape[0] + 1
            )
            return statistic, pvalue

        def __gamma_fit_pvalue(core):
            chrom, value, n_bins = core["chromosome"], core["value"], core["n_bins"]
            statistic = np.abs(value) * n_bins
            shape, loc, scale = self._null_dist["gamma-fit"][chrom]
            pvalue = 1 - gamma.cdf(statistic, shape, loc, scale)
            return statistic, pvalue

        def __delta_pvalue(core):
            chrom, v = core["chromosome"], core["binned_delta_ranks"]
            if (~np.isnan(v)).sum() > 0:
                statistic = np.max(np.abs(v[~np.isnan(v)]))
                pvalue = ((self._null_dist["delta"][chrom] > statistic).sum() + 1) / (
                    self._null_dist["delta"][chrom].shape[0] + 1
                )
            else:
                statistic = np.NaN
                pvalue = np.NaN
            return statistic, pvalue

        COMPUTE_PVALUE = {
            "gamma": __gamma_pvalue,
            "gamma-fit": __gamma_fit_pvalue,
            "delta": __delta_pvalue,
        }

        if self._null_dist is None:
            return None
        else:
            result = {}
            for test in self._stat_test:
                statistic, pvalue = COMPUTE_PVALUE[test](core)
                result[test] = (statistic, pvalue)
            return result

    def segment_chromosome(
        self, s1: CalderSubCompartments, s2: CalderSubCompartments
    ) -> CalderDifferentialCompartments:
        s12_rank = (
            self.get_binned_deltaRank(s1, s2)
            .sort_values(["chr", "start", "end"])
            .reset_index(drop=True)
        )
        values = s12_rank["delta_rank"].values

        expected_std = {1: np.inf, values.shape[0]: 0}

        def __get_expected_std(v, n_points=100):
            if v.shape[0] in expected_std.keys():
                return expected_std[v.shape[0]]
            else:
                stds = np.zeros(n_points, dtype=float)
                for p in range(n_points):
                    start_p = self._random_generator.choice(
                        values.shape[0] - v.shape[0]
                    )
                    end_p = start_p + v.shape[0]
                    seg = values[start_p:end_p]
                    if np.sum(~np.isnan(seg)) > 0:
                        std_p = np.nanstd(seg)
                        stds[p] = std_p
                expected_std[v.shape[0]] = np.nanmean(stds)
                return expected_std[v.shape[0]]

        def __pos_neg_segmentation(v, segment_value_function=np.mean):
            """Given a signal, it segments it based on positive and negative values
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
                seg_values.append(segment_value_function(v[seg_start : i + 1]))
            segments = np.array(segments, dtype=int)
            seg_values = np.array(seg_values)
            return segments, seg_values

        def __recursive_segmentation(X, sub_mean=True):
            v = X["delta_rank"].values
            rank1 = np.nanmean(X["domain_rank_1"].values)
            rank2 = np.nanmean(X["domain_rank_2"].values)
            mean_signal = np.nanmean(v)
            std_signal = np.nanstd(v)
            exp_std = __get_expected_std(v) if self._min_std is None else self._min_std
            if (std_signal >= exp_std) or (np.isnan(v).sum() > 0):
                segments, _ = __pos_neg_segmentation(v - mean_signal if sub_mean else v)
                all_rrs = []
                for s in range(segments.shape[0]):
                    seg = segments[s, :]
                    Xr = X.iloc[seg[0] : seg[1]]
                    rrs = __recursive_segmentation(Xr)
                    all_rrs.append(rrs)
                result = pd.concat(all_rrs, axis=0, ignore_index=True)
            else:
                if X.delta_rank.hasnans:
                    result = bioframe.merge(
                        X.dropna(
                            subset=["delta_rank"].rename(columns={"chr": "chrom"})
                        )[["chrom", "start", "end"]]
                    ).rename(columns={"chrom": "chr"})
                else:
                    result = pd.DataFrame.from_dict(
                        [
                            {
                                "chr": X.chr.iloc[0],
                                "start": X.start.min(),
                                "end": X.end.max(),
                            }
                        ]
                    )
                result["value"] = mean_signal
                result["rank1"] = rank1
                result["comp1"] = s1.discretize_rank(rank1, chrom=result.chr.iloc[0])
                result["rank2"] = rank2
                result["comp2"] = s2.discretize_rank(rank2, chrom=result.chr.iloc[0])
                stat_result = self.get_pvalue(
                    core={
                        "chromosome": X.chr.iloc[0],
                        "binned_delta_ranks": v,
                        "value": mean_signal,
                        "length": result.end.max() - result.start.min(),
                        "n_bins": (result.end.max() - result.start.min())
                        // self.binSize,
                    }
                )
                if stat_result is not None:
                    for test in stat_result.keys():
                        result[f"{test}_statistic"] = stat_result[test][0]
                        result[f"{test}_pvalue"] = stat_result[test][1]
            return result

        result = __recursive_segmentation(s12_rank, sub_mean=False)
        expected_std = [
            {
                "chr": s12_rank.chr.iloc[0],
                "segment_size": segment_size,
                "expected_std": exp_std,
            }
            for segment_size, exp_std in expected_std.items()
        ]
        expected_std = pd.DataFrame.from_dict(expected_std)
        expected_std = expected_std.sort_values(["chr", "segment_size"]).reset_index(
            drop=True
        )
        return CalderDifferentialCompartments(
            input_df=result,
            signal_df=s12_rank[
                [
                    "chr",
                    "start",
                    "end",
                    "delta_rank",
                    "domain_rank_1",
                    "domain_rank_2",
                    "compartment_label_8_1",
                    "compartment_label_8_2",
                ]
            ],
            expected_std=expected_std,
        )

