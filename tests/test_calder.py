import pytest
from diffcomp.calder import CalderSubCompartments, COMPS_HEADER

__author__ = "lucananni93"
__copyright__ = "lucananni93"
__license__ = "MIT"

FILE1="tests/data/calder/RPE_TP53_Ctrl_mega_CompartmentDomains.bed"


def test_read():
	comps = CalderSubCompartments.read(FILE1)
	assert comps.columns.tolist() == COMPS_HEADER

def test_domains():
	comps = CalderSubCompartments(FILE1)
	assert comps.domains.columns.tolist() == COMPS_HEADER

def test_binnify():
	comps = CalderSubCompartments(FILE1)
	binned_comps = comps.binnify(50000)
	assert binned_comps.iloc[0].end - binned_comps.iloc[0].start == 50000

	binned_comps2 = comps.binnify(10000)
	assert binned_comps2.iloc[0].end - binned_comps2.iloc[0].start == 10000

	assert sorted(comps._binned.keys()) == [10000, 50000]

def test_get_chromosome():
	comps = CalderSubCompartments(FILE1)
	comps_chr = comps.get_chromosomes("chr1")
	assert list(comps_chr.domains.chr.unique()) == ['chr1']
