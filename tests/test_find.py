import pytest
import pandas as pd
import tempfile
import os
from unittest.mock import patch
from io import StringIO

from vmwhere.find import (
    find_motif_genomic_coordinates,
    merge_adjacent_motifs,
    name_each_region,
    run_find
)


class TestFindMotifGenomicCoordinates:
    """Tests for find_motif_genomic_coordinates function."""

    def test_single_motif_occurrence(self):
        """Test finding a single occurrence of a motif."""
        # Sequence: A T C G G G A A A T C G
        # Index:    0 1 2 3 4 5 6 7 8 9 10 11
        # GGAA starts at position 4 (G G A A at 4,5,6,7)
        sequence = "ATCGGGAAATCG"
        motif = "GGAA"
        buffer_size = 2

        results = find_motif_genomic_coordinates(sequence, motif, buffer_size)

        assert len(results) == 1
        assert results[0]['Motif'] == "GGAA"
        assert results[0]['Motif_Start'] == 4
        assert results[0]['Motif_End'] == 8
        assert results[0]['Perfect_Repeats'] == 1
        assert results[0]['Total_Repeats'] == 1

    def test_multiple_consecutive_repeats(self):
        """Test finding multiple consecutive repeats of a motif."""
        sequence = "ATCGGAAGGAAGGAAATCG"
        motif = "GGAA"
        buffer_size = 2

        results = find_motif_genomic_coordinates(sequence, motif, buffer_size)

        assert len(results) == 1
        assert results[0]['Perfect_Repeats'] == 3
        assert results[0]['Total_Repeats'] == 3
        assert results[0]['Motif_Start'] == 3
        assert results[0]['Motif_End'] == 15

    def test_multiple_separate_occurrences(self):
        """Test finding multiple separate occurrences of a motif."""
        sequence = "GGAAATCGATCGGGAA"
        motif = "GGAA"
        buffer_size = 2

        results = find_motif_genomic_coordinates(sequence, motif, buffer_size)

        assert len(results) == 2
        assert results[0]['Motif_Start'] == 0
        assert results[1]['Motif_Start'] == 12

    def test_no_motif_found(self):
        """Test when motif is not present in sequence."""
        sequence = "ATCGATCGATCG"
        motif = "GGAA"
        buffer_size = 2

        results = find_motif_genomic_coordinates(sequence, motif, buffer_size)

        assert len(results) == 0

    def test_buffer_at_sequence_start(self):
        """Test buffer calculation when motif is at sequence start."""
        sequence = "GGAAGGAAATCG"
        motif = "GGAA"
        buffer_size = 10

        results = find_motif_genomic_coordinates(sequence, motif, buffer_size)

        assert len(results) == 1
        # Buffer start should be clamped to 0
        assert results[0]['Buffer_Start'] == 0
        assert results[0]['Buffer_End'] == 18  # 8 + 10

    def test_buffer_calculation(self):
        """Test buffer calculation around motif."""
        sequence = "ATCGATCGGGAAGGAAATCGATCG"
        motif = "GGAA"
        buffer_size = 3

        results = find_motif_genomic_coordinates(sequence, motif, buffer_size)

        assert len(results) == 1
        assert results[0]['Buffer_Start'] == 5  # 8 - 3
        assert results[0]['Buffer_End'] == 19  # 16 + 3

    def test_empty_sequence(self):
        """Test with empty sequence."""
        sequence = ""
        motif = "GGAA"
        buffer_size = 2

        results = find_motif_genomic_coordinates(sequence, motif, buffer_size)

        assert len(results) == 0

    def test_motif_at_end(self):
        """Test when motif is at the end of sequence."""
        sequence = "ATCGATCGGGAA"
        motif = "GGAA"
        buffer_size = 2

        results = find_motif_genomic_coordinates(sequence, motif, buffer_size)

        assert len(results) == 1
        assert results[0]['Motif_End'] == 12

    def test_overlapping_motif_regions(self):
        """Test behavior with potentially overlapping motifs (non-overlapping search)."""
        # GGAAGGAA has 2 GGAA motifs that are consecutive
        sequence = "ATCGGAAGGAAATCG"
        motif = "GGAA"
        buffer_size = 0

        results = find_motif_genomic_coordinates(sequence, motif, buffer_size)

        # Should find them as one consecutive region
        assert len(results) == 1
        assert results[0]['Perfect_Repeats'] == 2


class TestMergeAdjacentMotifs:
    """Tests for merge_adjacent_motifs function."""

    def test_merge_close_regions(self):
        """Test merging regions that are within max_gap."""
        df = pd.DataFrame({
            'chrom': ['chr1', 'chr1'],
            'motif': ['GGAA', 'GGAA'],
            'Motif_Start': [100, 120],
            'Motif_End': [110, 130],
            'Buffer_Start': [90, 110],
            'Buffer_End': [120, 140],
            'Total_Repeats': [2, 3],
            'Perfect_Repeats': [2, 3]
        })

        merged = merge_adjacent_motifs(df, max_gap=15)

        assert len(merged) == 1
        assert merged.iloc[0]['Motif_Start'] == 100
        assert merged.iloc[0]['Motif_End'] == 130
        assert merged.iloc[0]['Total_Repeats'] == 5
        assert merged.iloc[0]['Perfect_Repeats'] == 5

    def test_no_merge_distant_regions(self):
        """Test that distant regions are not merged."""
        df = pd.DataFrame({
            'chrom': ['chr1', 'chr1'],
            'motif': ['GGAA', 'GGAA'],
            'Motif_Start': [100, 200],
            'Motif_End': [110, 210],
            'Buffer_Start': [90, 190],
            'Buffer_End': [120, 220],
            'Total_Repeats': [2, 3],
            'Perfect_Repeats': [2, 3]
        })

        merged = merge_adjacent_motifs(df, max_gap=15)

        assert len(merged) == 2

    def test_no_merge_different_chromosomes(self):
        """Test that regions on different chromosomes are not merged."""
        df = pd.DataFrame({
            'chrom': ['chr1', 'chr2'],
            'motif': ['GGAA', 'GGAA'],
            'Motif_Start': [100, 105],
            'Motif_End': [110, 115],
            'Buffer_Start': [90, 95],
            'Buffer_End': [120, 125],
            'Total_Repeats': [2, 3],
            'Perfect_Repeats': [2, 3]
        })

        merged = merge_adjacent_motifs(df, max_gap=50)

        # Should NOT merge even though positions are close, because different chromosomes
        assert len(merged) == 2
        assert set(merged['chrom'].unique()) == {'chr1', 'chr2'}

    def test_empty_dataframe(self):
        """Test with empty dataframe."""
        df = pd.DataFrame(columns=['chrom', 'motif', 'Motif_Start', 'Motif_End',
                                   'Buffer_Start', 'Buffer_End', 'Total_Repeats', 'Perfect_Repeats'])

        merged = merge_adjacent_motifs(df, max_gap=10)

        assert len(merged) == 0

    def test_single_region(self):
        """Test with single region (nothing to merge)."""
        df = pd.DataFrame({
            'chrom': ['chr1'],
            'motif': ['GGAA'],
            'Motif_Start': [100],
            'Motif_End': [110],
            'Buffer_Start': [90],
            'Buffer_End': [120],
            'Total_Repeats': [2],
            'Perfect_Repeats': [2]
        })

        merged = merge_adjacent_motifs(df, max_gap=10)

        assert len(merged) == 1

    def test_merge_multiple_regions(self):
        """Test merging multiple consecutive regions."""
        df = pd.DataFrame({
            'chrom': ['chr1', 'chr1', 'chr1'],
            'motif': ['GGAA', 'GGAA', 'GGAA'],
            'Motif_Start': [100, 115, 130],
            'Motif_End': [110, 125, 140],
            'Buffer_Start': [90, 105, 120],
            'Buffer_End': [120, 135, 150],
            'Total_Repeats': [2, 2, 2],
            'Perfect_Repeats': [2, 2, 2]
        })

        merged = merge_adjacent_motifs(df, max_gap=10)

        assert len(merged) == 1
        assert merged.iloc[0]['Total_Repeats'] == 6
        assert merged.iloc[0]['Perfect_Repeats'] == 6

    def test_buffer_bounds_after_merge(self):
        """Test that buffer bounds are correctly calculated after merge."""
        df = pd.DataFrame({
            'chrom': ['chr1', 'chr1'],
            'motif': ['GGAA', 'GGAA'],
            'Motif_Start': [100, 120],
            'Motif_End': [110, 130],
            'Buffer_Start': [80, 110],
            'Buffer_End': [120, 150],
            'Total_Repeats': [2, 3],
            'Perfect_Repeats': [2, 3]
        })

        merged = merge_adjacent_motifs(df, max_gap=15)

        assert len(merged) == 1
        # Buffer_Start should be min of all buffer starts
        assert merged.iloc[0]['Buffer_Start'] == 80
        # Buffer_End should be max of all buffer ends
        assert merged.iloc[0]['Buffer_End'] == 150


class TestNameEachRegion:
    """Tests for name_each_region function."""

    def test_basic_naming(self):
        """Test basic region naming."""
        df = pd.DataFrame({
            'chrom': ['chr1', 'chr1', 'chr2'],
            'Buffer_Start': [100, 200, 100],
            'Buffer_End': [150, 250, 150],
            'motif': ['GGAA', 'GGAA', 'GGAA']
        })

        named = name_each_region(df)

        assert 'region_name' in named.columns
        assert named.iloc[0]['region_name'] == 'chr1_region_1'
        assert named.iloc[1]['region_name'] == 'chr1_region_2'
        assert named.iloc[2]['region_name'] == 'chr2_region_3'

    def test_naming_preserves_original(self):
        """Test that naming doesn't modify original dataframe."""
        df = pd.DataFrame({
            'chrom': ['chr1'],
            'Buffer_Start': [100],
            'Buffer_End': [150],
            'motif': ['GGAA']
        })

        original_columns = list(df.columns)
        name_each_region(df)

        assert list(df.columns) == original_columns

    def test_empty_dataframe(self):
        """Test with empty dataframe."""
        df = pd.DataFrame(columns=['chrom', 'Buffer_Start', 'Buffer_End', 'motif'])

        named = name_each_region(df)

        assert 'region_name' in named.columns
        assert len(named) == 0


class TestRunFind:
    """Integration tests for run_find function."""

    def test_basic_run_find(self):
        """Test basic run_find with simple FASTA file."""
        # Create a temporary FASTA file
        fasta_content = ">chr1\nATCGGGAAGGAAGGAAATCGATCGGGAAGGAAATCG\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, "test.fasta")
            with open(fasta_path, 'w') as f:
                f.write(fasta_content)

            output_dir = os.path.join(tmpdir, "output")

            result = run_find(
                motif="GGAA",
                repeats=2,
                gap=10,
                fasta_file=fasta_path,
                buffer=5,
                output_dir=output_dir
            )

            # Check output file was created
            output_file = os.path.join(output_dir, "microsatellite_coordinates.bed")
            assert os.path.exists(output_file)

            # Check result dataframe
            assert result is not None
            assert len(result) > 0

    def test_run_find_no_matches(self):
        """Test run_find when no motifs are found."""
        fasta_content = ">chr1\nATCGATCGATCGATCG\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, "test.fasta")
            with open(fasta_path, 'w') as f:
                f.write(fasta_content)

            output_dir = os.path.join(tmpdir, "output")

            result = run_find(
                motif="GGAA",
                repeats=2,
                gap=10,
                fasta_file=fasta_path,
                buffer=5,
                output_dir=output_dir
            )

            # Should return None when no motifs found
            assert result is None

    def test_run_find_reverse_complement(self):
        """Test that reverse complement motifs are also found."""
        # GGAA reverse complement is TTCC
        fasta_content = ">chr1\nATCGTTCCTTCCTTCCATCG\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, "test.fasta")
            with open(fasta_path, 'w') as f:
                f.write(fasta_content)

            output_dir = os.path.join(tmpdir, "output")

            result = run_find(
                motif="GGAA",
                repeats=2,
                gap=10,
                fasta_file=fasta_path,
                buffer=5,
                output_dir=output_dir
            )

            # Should find the reverse complement
            assert result is not None
            assert len(result) > 0

    def test_run_find_multiple_chromosomes(self):
        """Test run_find with multiple chromosomes."""
        fasta_content = ">chr1\nATCGGGAAGGAAGGAAATCG\n>chr2\nGGAAGGAAGGAAGGAA\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, "test.fasta")
            with open(fasta_path, 'w') as f:
                f.write(fasta_content)

            output_dir = os.path.join(tmpdir, "output")

            result = run_find(
                motif="GGAA",
                repeats=2,
                gap=10,
                fasta_file=fasta_path,
                buffer=5,
                output_dir=output_dir
            )

            assert result is not None
            # Should have regions from both chromosomes
            assert len(result['chrom'].unique()) == 2

    def test_run_find_filters_by_repeats(self):
        """Test that run_find filters by minimum perfect repeats."""
        # One region with 2 repeats, one with 4
        fasta_content = ">chr1\nGGAAGGAAATCGATCGATCGATCGATCGGGAAGGAAGGAAGGAA\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, "test.fasta")
            with open(fasta_path, 'w') as f:
                f.write(fasta_content)

            output_dir = os.path.join(tmpdir, "output")

            # Require at least 3 repeats
            result = run_find(
                motif="GGAA",
                repeats=3,
                gap=5,
                fasta_file=fasta_path,
                buffer=5,
                output_dir=output_dir
            )

            # Should only find the region with 4 repeats
            assert result is not None
            assert len(result) == 1

    def test_run_find_output_format(self):
        """Test that output file has correct BED format."""
        fasta_content = ">chr1\nATCGGGAAGGAAGGAAATCG\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, "test.fasta")
            with open(fasta_path, 'w') as f:
                f.write(fasta_content)

            output_dir = os.path.join(tmpdir, "output")

            run_find(
                motif="GGAA",
                repeats=2,
                gap=10,
                fasta_file=fasta_path,
                buffer=5,
                output_dir=output_dir
            )

            # Read and verify output format
            output_file = os.path.join(output_dir, "microsatellite_coordinates.bed")
            with open(output_file, 'r') as f:
                lines = f.readlines()

            assert len(lines) > 0
            # BED format should have tab-separated columns
            fields = lines[0].strip().split('\t')
            assert len(fields) == 5  # chrom, start, end, motif, region_name

    def test_run_find_case_insensitive_motif(self):
        """Test that motif search is case-insensitive."""
        fasta_content = ">chr1\nATCGggaaggaaggaaATCG\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, "test.fasta")
            with open(fasta_path, 'w') as f:
                f.write(fasta_content)

            output_dir = os.path.join(tmpdir, "output")

            result = run_find(
                motif="ggaa",  # lowercase
                repeats=2,
                gap=10,
                fasta_file=fasta_path,
                buffer=5,
                output_dir=output_dir
            )

            assert result is not None
            assert len(result) > 0


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_very_long_motif(self):
        """Test with a longer motif."""
        sequence = "ATCGGGAATTCCGGAATTCCGGAATTCCATCG"
        motif = "GGAATTCC"
        buffer_size = 2

        results = find_motif_genomic_coordinates(sequence, motif, buffer_size)

        assert len(results) == 1
        assert results[0]['Perfect_Repeats'] == 3

    def test_single_nucleotide_motif(self):
        """Test with single nucleotide motif (homopolymer)."""
        # Use a sequence where the homopolymer is the only occurrence
        sequence = "TCGAAAAAAAATCG"
        motif = "A"
        buffer_size = 2

        results = find_motif_genomic_coordinates(sequence, motif, buffer_size)

        assert len(results) == 1
        assert results[0]['Perfect_Repeats'] == 8

    def test_motif_same_as_sequence(self):
        """Test when motif equals the entire sequence."""
        sequence = "GGAAGGAAGGAA"
        motif = "GGAA"
        buffer_size = 0

        results = find_motif_genomic_coordinates(sequence, motif, buffer_size)

        assert len(results) == 1
        assert results[0]['Perfect_Repeats'] == 3
        assert results[0]['Motif_Start'] == 0
        assert results[0]['Motif_End'] == 12

    def test_zero_buffer_size(self):
        """Test with zero buffer size."""
        sequence = "ATCGGGAAATCG"
        motif = "GGAA"
        buffer_size = 0

        results = find_motif_genomic_coordinates(sequence, motif, buffer_size)

        assert len(results) == 1
        assert results[0]['Buffer_Start'] == results[0]['Motif_Start']
        assert results[0]['Buffer_End'] == results[0]['Motif_End']

    def test_large_buffer_size(self):
        """Test with buffer size larger than sequence."""
        sequence = "GGAA"
        motif = "GGAA"
        buffer_size = 1000

        results = find_motif_genomic_coordinates(sequence, motif, buffer_size)

        assert len(results) == 1
        assert results[0]['Buffer_Start'] == 0  # Clamped to 0
        assert results[0]['Buffer_End'] == 1004  # 4 + 1000
