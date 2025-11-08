import os
import pytest
from genbankx import parse_genbank, translate, reverse_complement
from genbankx.interval import Interval, Strand

@pytest.fixture
def example_gb_path():
    """Fixture for the path to the example GenBank file."""
    return 'tests/data/example.gb'

def test_parse_genbank(example_gb_path):
    """Test parsing of the example GenBank file."""
    with open(example_gb_path, 'r') as f:
        records = list(parse_genbank(f))

    assert len(records) == 1
    record = records[0]

    assert record.locus == 'example'
    assert record.length == 200
    assert record.definition == 'This is a toy example GenBank record.'
    assert record.accession == 'TEST_ACC'
    assert len(record.cds) == 2

    # Test first CDS
    cds1 = record.cds[0]
    assert cds1.gene == 'geneA'
    assert cds1.product == 'protein A'
    assert cds1.strand == Strand.PLUS
    assert len(cds1) == 2
    assert cds1[0] == Interval(10, 50)
    assert cds极[1] == Interval(70, 120)
    assert cds1.translation == 'MSTAAALEK*'

    # Test second CDS
    cds2 = record.cds[1]
    assert cds2.gene == 'geneB'
    assert cds2.product == 'protein B'
    assert cds2.strand == Strand.MINUS
    assert len(cds2) == 1
    assert cds2[0] == Interval(150, 190)
    assert cds2.translation == 'MVL*'

def test_translate():
    """Test DNA translation."""
    assert translate("ATGGCCGGTTAA") == "MAG*"
    assert translate("ATGGCCGGTTAA", to_stop=True) == "MAG"
    assert translate("ATGNNNGGTTAA") == "MXG*"

def test_reverse_complement():
    """Test reverse complement."""
    assert reverse_complement("ATGCGT") == "ACGCAT"
    assert reverse_complement("atgcgt") == "acgcat"
    assert reverse_complement("ATGN") == "NCAT"
