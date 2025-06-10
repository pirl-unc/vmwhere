import re
import pytest
from src.pipeline import (
    decompose_seq_with_motif_anchors,
    recompose_string_from_structure
)


@pytest.mark.parametrize("seq, motif, expected_prefixes", [
    ("GGAA" * 4, "GGAA", ["4GGAA"]),
    ("GGAA" * 3 + "TC" + "GGAA" * 2, "GGAA", ["3GGAA", "1TC", "2GGAA"]),
    ("GGAA" * 2 + "AGAA" + "GGAA", "GGAA", ["2GGAA", "1AGAA", "1GGAA"]),
    ("TTCGGAA" * 2 + "GGAC" + "GGAA" * 2, "GGAA", ["1GGAA", "1TTC", "1GGAA" , "1GGAC", "2GGAA"]),
    ("TCGGGATGGAAGGAA", "GGAA", ["2GGAA"]),
    ("GGGAGGGAGGGAAGGAAGGAGGGAAGGGAAGAAGGAAGGAAGGAG" + "GGAA" * 12 + "GGAG", "GGAA", ["2GGGA", "1G", "2GGAA", "1GGAG", "1GGAA", "1G", "1GGAA", "1GAA", "2GGAA", "1GGAG", "12GGAA"]),
    ("GGATGGATGGATGGAT", "GGAA", [None]),  # no motif present
])

def test_decompose_and_recompose(seq, motif, expected_substructure):
    decomp_structure, length, subseq = decompose_seq_with_motif_anchors(seq, motif)

    # structure might be None if motif isn't found
    if structure is None:
        assert expected_prefixes == [None]
    else:
        # Check expected structure components
        for decomp_substructure in expected_substructure:
            assert prefix in structure

        # Check that the decomposition didn't add or remove any components by reconstructing seq from the decomposed sequence and comparing that to the original sequence
        recomposed = recompose_string_from_structure(structure)
        assert recomposed == subseq

