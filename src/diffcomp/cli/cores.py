import argparse
import logging
import sys
import os

from .. import __version__
from ..calder import CalderSubCompartments
from ..diff import CalderRecursiveDifferentialSegmentator


__author__ = "Luca Nanni"
__contact__ = "luca.nanni@unil.ch"
__date__ = "2022/02/24"

_logger = logging.getLogger(__name__)



def parse_args():
    parser = argparse.ArgumentParser(description="""
        Identifying Compartment Repositioning Events from Calder genomic segmentations.
        Given sample1 and sample2 Calder segmentations, it identifies regions undergoing statistically
        significant compartment repositioning in sample2 in comparison to sample1.

        Significance of the repositioning is determined using paired control samples (control2 vs control1),
        which are provided by the user. Usually, replicates of the same experiments are used to model the
        intrinsic biological noise in the Hi-C compartment calls.
        """)
    parser.add_argument("sample1_path", type=str, help="Path to the Calder segmentation of sample 1")
    parser.add_argument("sample2_path", type=str, help="Path to the Calder segmentation of sample 2")
    parser.add_argument("binsize", type=int, help="Resolution to use in the analysis")
    parser.add_argument("output_path", type=str, help="Path where to store the identified regions")
    parser.add_argument("--algo", type=str, help="Which algorithm to use for finding CoREs", default='recursive')
    parser.add_argument("--statistical_test", type=str, help="Which test to use to determine significance of CoREs (all, gamma, gamma-fit, delta)", default='gamma')
    parser.add_argument("--min_std", type=float, default=None, help="Maximum standard deviation allowed for segmented regions (DEPRECATED)")
    parser.add_argument("--control1_path", type=str, help="Path(s) to the Calder segmentation(s) to use to use as control 1 (comma-separated)")
    parser.add_argument("--control2_path", type=str, help="Path(s) to the Calder segmentation(s) to use to use as control 2 (comma-separated)")
    parser.add_argument("--signal_path", type=str, help = "Path where to store the binned differential signal")
    parser.add_argument("--bed_path", type=str, help = "Path where to store the identified regions in BED format")
    parser.add_argument("--control_distribution_path", type=str, help = "Path where to store the computed control distribution (Pickle format)")
    parser.add_argument("--coordinates", type=str, default="zero-based", help="Coordinate system of the input files (zero-based / one-based)")
    parser.add_argument("--genome", type=str, default="hg19", help="Genome (Default: hg19)")
    parser.add_argument("--chromosomes", type=str, default="ALL", help="List of chromosomes to perform the analysis on (Default: all chromosomes, comma-separated)")
    parser.add_argument("--verbose", dest="loglevel", help="Set loglevel to INFO", action="store_const", const=logging.INFO)
    parser.add_argument("--very-verbose", dest="loglevel", help="Set loglevel to DEBUG", action="store_const", const=logging.DEBUG)
    parser.add_argument("--version", action="version",version="diffComp {ver}".format(ver=__version__))
    return parser.parse_args()


def setup_logging(loglevel):
    logformat = "[%(asctime)s] %(levelname)s:%(name)s:%(message)s"
    logging.basicConfig(level=loglevel, stream=sys.stdout, format=logformat, datefmt="%Y-%m-%d %H:%M:%S")


def main():
    args = parse_args()
    setup_logging(args.loglevel)

    if args.algo == 'recursive':
        _logger.info("Recursive segmentation algorithm was selected")
        _logger.info(f"- binSize = {args.binsize}")
        _logger.info(f"- min_std = {args.min_std}")
        _logger.info(f"- statistical test = {args.statistical_test}")
        _logger.info(f"- Control distribution path = {args.control_distribution_path}")
        segmentator = CalderRecursiveDifferentialSegmentator(binSize = args.binsize, min_std = args.min_std, statistical_test = args.statistical_test)
    else:
        raise ValueError("Unknown CoREs algorithm")

    _logger.debug(f"Loading Calder compartments from {args.sample1_path}")
    comps1 = CalderSubCompartments(args.sample1_path, genome=args.genome, coordinates = args.coordinates)
    _logger.debug(f"Loading Calder compartments from {args.sample2_path}")
    comps2 = CalderSubCompartments(args.sample2_path, genome=args.genome, coordinates = args.coordinates)

    if args.chromosomes != "ALL":
        comps1 = comps1.get_chromosomes(args.chromosomes.split(","))
        comps2 = comps2.get_chromosomes(args.chromosomes.split(","))

    if os.path.exists(args.control_distribution_path):
        _logger.info("Loading control distribution from file")
        segmentator.load_control_distribution(args.control_distribution_path)
    else:
        if (args.control1_path is not None) and (args.control2_path is not None):
            args.control1_path = args.control1_path.split(",")
            args.control2_path = args.control2_path.split(",")
            if len(args.control1_path) != len(args.control2_path):
                raise ValueError("Control samples must be the same number in 'control1_path' and 'control2_path'")
            if args.algo != 'recursive':
                raise ValueError("Control samples are required only for the recursive segmentation algorithm")
            _logger.info("Building null distribution from control samples")
            controls = [(CalderSubCompartments(x, genome=args.genome, coordinates = args.coordinates),
                         CalderSubCompartments(y, genome=args.genome, coordinates = args.coordinates)) \
                            for x, y in zip(args.control1_path, args.control2_path)]
            if args.chromosomes != "ALL":
                controls = [(t[0].get_chromosomes(args.chromosomes.split(",")), t[1].get_chromosomes(args.chromosomes.split(","))) for t in controls]
            segmentator.build_control_distribution(controls)
            if args.control_distribution_path is not None:
                segmentator.write_control_distribution(args.control_distribution_path)

    _logger.info("Starting the segmentation")
    segments = segmentator.segment(comps1, comps2)

    _logger.info(f"Saving result to {args.output_path}")
    segments.to_tsv(args.output_path)

    if args.signal_path is not None:
        segments.write_signal(args.signal_path)
    if args.bed_path is not None:
        segments.to_bed(args.bed_path)


if __name__ == '__main__':
    main()
