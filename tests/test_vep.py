import pytest

from .test_vcf import valid_vcf, valid_vcf_file
from veppy.vep import VEP


@pytest.fixture
def vep(valid_vcf):
    vep_description = valid_vcf.info["CSQ"].Description
    return VEP(vep_description)


def test_vep(vep):
    assert all(
        x in vep.variables for x in ("pop1_af", "pop2_af", "SYMBOL", "Consequence")
    )


def test_vep_variant(valid_vcf, vep):
    v = next(valid_vcf.variants)
    vv = vep.variant(v)
    assert all(x in vv for x in ("pop1_af", "pop2_af", "SYMBOL", "Consequence"))
    assert vv["pop1_af"] == "0.00345"
    assert vv["pop2_af"] is None
    assert vv["SYMBOL"] == "YFG"


def test_vcf_annotations(valid_vcf):
    vep = valid_vcf.annotations
    v = next(valid_vcf.variants)
    vv = vep.variant(v)
    assert vv["pop1_af"] == "0.00345"
    assert vv["pop2_af"] is None
    assert vv["SYMBOL"] == "YFG"


def test_variant_annotation(valid_vcf):
    v = next(valid_vcf.variants)
    assert v.annotation["pop1_af"] == "0.00345"
    assert v.annotation["pop2_af"] is None
    assert v.annotation["SYMBOL"] == "YFG"


def test_missing_annotation(valid_vcf):
    next(valid_vcf.variants)
    v = next(valid_vcf.variants)
    assert v.annotation is None


def test_variable_length(valid_vcf):
    assert len(valid_vcf.annotations.variables) == 4


def test_variable_quoting(valid_vcf):
    assert valid_vcf.annotations.variables[-1][-1] != '"'
