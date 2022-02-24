import argparse
import logging
import sys

from .. import __version__
from ..calder import CalderSubCompartments
from ..diff import CalderRecursiveDifferentialSegmentator


__author__ = "Luca Nanni"
__contact__ = "luca.nanni@unil.ch"
__date__ = "2022/02/24"

_logger = logging.getLogger(__name__)



def parse_args():
    parser = argparse.ArgumentParser(description="Identifying Compartment Repositioning Events from Calder genomic segmentations")
    parser.add_argument(
        "--version",
        action="version",
        version="diffComp {ver}".format(ver=__version__),
    )
    parser.add_argument("sample1_path", type=str, help="Path to the Calder segmentation of sample 1")
    parser.add_argument("sample2_path", type=str, help="Path to the Calder segmentation of sample 2")
    parser.add_argument("binsize", type=int, help="Resolution to use in the analysis")
    parser.add_argument("output_path", type=str, help="Path where to store the identified regions")
    parser.add_argument("--algo", type=str, help="Which algorithm to use for finding CoREs", default='recursive')
    parser.add_argument("--min_std", type=float, default=0.1, help="Maximum standard deviation allowed for segmented regions")
    parser.add_argument("--control1_path", type=str, nargs='*', help="Path(s) to the Calder segmentation(s) to use to use as control 1")
    parser.add_argument("--control2_path", type=str, nargs='*', help="Path(s) to the Calder segmentation(s) to use to use as control 2")
    parser.add_argument(
        "-v",
        "--verbose",
        dest="loglevel",
        help="set loglevel to INFO",
        action="store_const",
        const=logging.INFO,
    )
    parser.add_argument(
        "-vv",
        "--very-verbose",
        dest="loglevel",
        help="set loglevel to DEBUG",
        action="store_const",
        const=logging.DEBUG,
    )
    return parser.parse_args()


def setup_logging(loglevel):
    logformat = "[%(asctime)s] %(levelname)s:%(name)s:%(message)s"
    logging.basicConfig(
        level=loglevel, stream=sys.stdout, format=logformat, datefmt="%Y-%m-%d %H:%M:%S"
    )


def main():
    args = parse_args()
    setup_logging(args.loglevel)

    if args.algo == 'recursive':
        _logger.info("Recursive segmentation algorithm was selected")
        _logger.info(f"- binSize = {args.binsize}")
        _logger.info(f"- min_std = {args.min_std}")
        segmentator = CalderRecursiveDifferentialSegmentator(binSize = args.binsize, min_std = args.min_std)
    else:
        raise ValueError("Unknown CoREs algorithm")

    _logger.debug(f"Loading Calder compartments from {args.sample1_path}")
    comps1 = CalderSubCompartments(args.sample1_path)
    _logger.debug(f"Loading Calder compartments from {args.sample2_path}")
    comps2 = CalderSubCompartments(args.sample2_path)


    if (args.control1_path is not None) and (args.control2_path is not None):
        if len(args.control1_path) != len(args.control2_path):
            raise ValueError("Control samples must be the same number in 'control1_path' and 'control2_path'")
        if args.algo != 'recursive':
            raise ValueError("Control samples are required only for the recursive segmentation algorithm")
        _logger.info("Building null distribution from control samples")
        controls = [(CalderSubCompartments(x), CalderSubCompartments(y)) for x, y in zip(args.control1_path, args.control2_path)]
        segmentator.build_control_distribution(controls)

    _logger.info("Starting the segmentation")
    segments = segmentator.segment(comps1, comps2)

    _logger.info(f"Saving result to {args.output_path}")
    segments.to_tsv(args.output_path)


if __name__ == '__main__':
    main()