import pytest
from diffcomp.calder import CalderSubCompartments
from diffcomp.diff import CalderDifferentialCompartments, \
						  CalderRecursiveDifferentialSegmentator

__author__ = "lucananni93"
__copyright__ = "lucananni93"
__license__ = "MIT"

CALDER_DIFF_COMP_FILE="tests/data/diff/RPE_TP53_Ctrl_mega_vs_RPE_TP53_20wT1_mega.tsv"
SEGMENTS_TSV_PATH="tests/data/diff/RPE_TP53_Ctrl_mega_vs_RPE_TP53_20wT1_chr20.tsv"
SEGMENTS_BED_PATH="tests/data/diff/RPE_TP53_Ctrl_mega_vs_RPE_TP53_20wT1_chr20.bed"
SIGNAL_PATH = "tests/data/diff/RPE_TP53_Ctrl_mega_vs_RPE_TP53_20wT1_chr20_signal.tsv"
FILE1="tests/data/calder/RPE_TP53_Ctrl_mega_CompartmentDomains.bed"
FILE2="tests/data/calder/RPE_TP53_20w0T1_mega_CompartmentDomains.bed"

CONTROL1="tests/data/calder/RPE_TP53_Ctrl_1_CompartmentDomains.bed"
CONTROL2="tests/data/calder/RPE_TP53_Ctrl_2_CompartmentDomains.bed"

def test_CalderDifferentialCompartments():
	seg = CalderDifferentialCompartments(path = CALDER_DIFF_COMP_FILE)

def test_segmentation():
	seg = CalderDifferentialCompartments(path = CALDER_DIFF_COMP_FILE)
	assert seg.segmentation.columns[: len(CalderDifferentialCompartments.CALDER_DIFFERENTIAL_COMPARTMENTS_HEADER)].tolist() == CalderDifferentialCompartments.CALDER_DIFFERENTIAL_COMPARTMENTS_HEADER


def test_recursive_segmentator():
	segmentator = CalderRecursiveDifferentialSegmentator(50000, random_state=42)
	assert segmentator.binSize == 50000
	comps1 = CalderSubCompartments(FILE1)
	comps2 = CalderSubCompartments(FILE2)
	comps12 = segmentator.get_binned_deltaRank(comps1, comps2)
	assert comps12.iloc[0].end - comps12.iloc[0].start == 50000

	control1 = CalderSubCompartments(CONTROL1)
	control2 = CalderSubCompartments(CONTROL2)

	null_dist = segmentator.build_control_distribution([(control1.get_chromosomes("chr20"), control2.get_chromosomes("chr20"))])
	assert set(null_dist["gamma"].keys()) == set(["chr20"])

	segments = segmentator.segment(comps1.get_chromosomes("chr20"), comps2.get_chromosomes("chr20"))
	assert segments.segmentation.columns[:len(CalderDifferentialCompartments.CALDER_DIFFERENTIAL_COMPARTMENTS_HEADER)].tolist() == CalderDifferentialCompartments.CALDER_DIFFERENTIAL_COMPARTMENTS_HEADER
	assert segments.signal.columns[:4].tolist() == ['chr', 'start', 'end', 'delta_rank']
	segments.to_tsv(SEGMENTS_TSV_PATH)
	segments.to_bed(SEGMENTS_BED_PATH)
	segments.write_signal(SIGNAL_PATH)
