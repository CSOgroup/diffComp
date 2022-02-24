import pytest
from diffcomp.calder import CalderSubCompartments
from diffcomp.diff import CalderDifferentialCompartments, CalderRecursiveDifferentialSegmentator

__author__ = "lucananni93"
__copyright__ = "lucananni93"
__license__ = "MIT"

CALDER_DIFF_COMP_FILE="tests/data/diff/TP53_20wTumo1_30_vs_RPE_Control_compartment_changes_recursive_segmentation.tsv"
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
	segmentator = CalderRecursiveDifferentialSegmentator(50000)
	assert segmentator.binSize == 50000
	comps1 = CalderSubCompartments(FILE1)
	comps2 = CalderSubCompartments(FILE2)
	comps12 = segmentator._merge_sample_compartments(comps1, comps2)
	assert comps12.iloc[0].end - comps12.iloc[0].start == 50000

	control1 = CalderSubCompartments(CONTROL1)
	control2 = CalderSubCompartments(CONTROL2)

	null_dist = segmentator.build_control_distribution([(control1, control2)])
	assert set(null_dist.keys()) == set(control1.domains.chr).intersection(set(control2.domains.chr))

	segments = segmentator.segment(comps1.get_chromosome("chr20"), comps2.get_chromosome("chr20"))
	assert segments.segmentation.columns.tolist() == CalderDifferentialCompartments.CALDER_DIFFERENTIAL_COMPARTMENTS_HEADER