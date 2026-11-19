import pandas as pd
from vmwhere.genotyper import identify_alleles, create_tsv_per_sample

REF = "GGAAGGAAGGAAGGAA"


def make_locus(region_id, alleles, supports):
    """One row per allele, in the order the genotyper would report them."""
    return pd.DataFrame({
        "CHROM": "chr1",
        "POS": 100,
        "ID": region_id,
        "REF": REF,
        "ALT": alleles,
        "END": 200,
        "MOTIF": "GGAA",
        "AL": [len(a) for a in alleles],
        "CN": [a.count("GGAA") for a in alleles],
        "CNM": [a.count("GGAA") for a in alleles],
        "MD": 1.0,
        "DS_READ": [f"{a.count('GGAA')}GGAA" for a in alleles],
        "DS_REF": "4GGAA",
        "RS": supports,
    })


def test_three_alleles_capped_to_two_best_supported():
    df = make_locus("r1", [REF, "GGAA" * 5, "GGAA" * 6], [12, 7, 11])
    called = identify_alleles(df, minor_thresh=0.2, homozygous_thresh=0.8)

    assert len(called) == 2
    # the two best-supported alleles, in their original order
    assert called["RS"].tolist() == [12, 11]
    assert called["ALT"].tolist() == [REF, "GGAA" * 6]


def test_tie_for_second_allele_keeps_first_listed():
    df = make_locus("r1", [REF, "GGAA" * 5, "GGAA" * 6], [12, 9, 9])
    called = identify_alleles(df, minor_thresh=0.2, homozygous_thresh=0.8)

    assert called["ALT"].tolist() == [REF, "GGAA" * 5]


def test_two_allele_and_homozygous_loci_unchanged():
    het = make_locus("het", ["GGAA" * 5, REF], [10, 14])
    hom = make_locus("hom", ["GGAA" * 6, REF], [19, 1])
    called = identify_alleles(pd.concat([het, hom], ignore_index=True),
                              minor_thresh=0.2, homozygous_thresh=0.8)

    het_called = called[called["ID"] == "het"]
    hom_called = called[called["ID"] == "hom"]
    assert het_called["ALT"].tolist() == ["GGAA" * 5, REF]
    assert hom_called["ALT"].tolist() == ["GGAA" * 6]


def test_capped_locus_reports_diploid_genotype():
    df = make_locus("r1", [REF, "GGAA" * 5, "GGAA" * 6, "GGAA" * 7], [8, 9, 7, 6])
    tsv = create_tsv_per_sample(identify_alleles(df, minor_thresh=0.2, homozygous_thresh=0.8))

    row = tsv.iloc[0]
    assert row["GT"] == "0/1"
    assert row["ALT"] == "GGAA" * 5
    for field in ["AL", "CN", "CNM", "MD", "DS_READ", "RS"]:
        assert len(row[field].split(",")) == 2
    assert row["RS"] == "8,9"
