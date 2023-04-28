import pytest

from veppy.vcf import VCF


@pytest.fixture
def valid_vcf_file(tmp_path):
    vcf_path = tmp_path / "valid.vcf"
    vcf_text = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=.,Type=String,Description="Sample genotype">
##FORMAT=<ID=AO,Number=.,Type=String,Description="Reads supporting each alternative allele, in order of appearance">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Alternative allele frequency">
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: SYMBOL|pop1_af|pop2_af|Consequence">
##INFO=<ID=CIGAR,Number=A,Type=String,Description="The extended CIGAR representation of each alternate allele, with the exception that '=' is replaced by 'M' to ease VCF parsing.  Note that INDEL alleles do not have the first matched base (which is provided by default, per the spec) referred to by the CIGAR.">
##INFO=<ID=DECOMPOSED,Number=.,Type=String,Description="The variant is a result of a decomposition of an indel or a complex variant">
##contig=<ID=chr1,length=100>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t10\t.\tA\tT\t100\tPASS\tCSQ=YFG|0.00345||123\tGT:AO:AF\t0/1:123:0.54
chr1\t15\t.\tC\tG\t123\tPASS\t\tGT:AO:AF\t0/1:123:0.54
chr1\t20\t.\tT\tG\t123\tPASS\tDECOMPOSED\tGT:AO:AF\t0/1:35:0.07
chr1\t20\t.\tT\tG\t.\tPASS\tDECOMPOSED\tGT:AO:AF\t0/1:35:0.07
"""
    with open(vcf_path, "w") as of:
        of.write(vcf_text)

    return vcf_path


@pytest.fixture
def valid_vcf(valid_vcf_file):
    return VCF(valid_vcf_file)


@pytest.fixture
def expected_info_keys():
    return [["CSQ"], [], ["DECOMPOSED"]]


@pytest.fixture
def expected_default_fields():
    return dict(
        CHROM=["chr1", "chr1", "chr1", "chr1"],
        POS=[10, 15, 20, 20],
        ID=[".", ".", ".", "."],
        REF=["A", "C", "T", "T"],
        ALT=[["T"], ["G"], ["G"], ["G"]],
        QUAL=[100.0, 123.0, 123.0, None],
        FILTER=["PASS", "PASS", "PASS", "PASS"],
    )


def test_vcf_init(valid_vcf_file):
    vcf = VCF(valid_vcf_file)
    assert vcf.samples == ["sample1"]
    assert vcf.version == "VCFv4.2"


def test_fields(valid_vcf, expected_default_fields):
    for i, v in enumerate(valid_vcf.variants):
        assert v.chrom == expected_default_fields["CHROM"][i]
        assert v["CHROM"] == expected_default_fields["CHROM"][i]
        assert v.pos == expected_default_fields["POS"][i]
        assert v["POS"] == expected_default_fields["POS"][i]
        assert v.id == expected_default_fields["ID"][i]
        assert v["ID"] == expected_default_fields["ID"][i]
        assert v.ref == expected_default_fields["REF"][i]
        assert v["REF"] == expected_default_fields["REF"][i]
        assert v.alt == expected_default_fields["ALT"][i]
        assert v["ALT"] == expected_default_fields["ALT"][i]
        assert v.qual == expected_default_fields["QUAL"][i]
        assert v["QUAL"] == expected_default_fields["QUAL"][i]
        assert v.filter == expected_default_fields["FILTER"][i]
        assert v["FILTER"] == expected_default_fields["FILTER"][i]


def test_format(valid_vcf):
    assert valid_vcf.format["AF"].ID == "AF"
    assert valid_vcf.format["AF"].Number == "A"
    assert valid_vcf.format["AF"].Type == "Float"
    assert valid_vcf.format["AF"].Description == "Alternative allele frequency"
    assert (
        valid_vcf.format["AO"].Description
        == "Reads supporting each alternative allele, in order of appearance"
    )


def test_info(valid_vcf):
    assert (
        valid_vcf.info["CIGAR"].Description
        == "The extended CIGAR representation of each alternate allele, with the exception that '=' is replaced by 'M' to ease VCF parsing.  Note that INDEL alleles do not have the first matched base (which is provided by default, per the spec) referred to by the CIGAR."
    )


def test_variant_info(valid_vcf, expected_info_keys):
    for v, ik in zip(valid_vcf.variants, expected_info_keys):
        for k in ik:
            assert k in v.info


def test_samples(valid_vcf):
    v = next(valid_vcf.variants)
    assert v.samples[0] == v.samples["sample1"]


def test_description_quoting(valid_vcf):
    assert valid_vcf.info["CSQ"].Description[-1] != '"'
