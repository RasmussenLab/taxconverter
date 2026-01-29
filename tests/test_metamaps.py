"""Tests for MetaMaps parser - validates user-observable behavior."""

import pytest
from taxconverter.metamaps import parse_metamaps_krona


class TestMetaMapsValidInput:
    """Test that valid MetaMaps input produces valid output."""

    def test_parse_valid_file(self, ncbi_ranks, fixtures_dir):
        """Parse valid metamaps file and validate all output properties."""
        path = fixtures_dir / "metamaps_valid.tsv"

        # Read input to count contigs (excluding comments and empty lines)
        with open(path) as f:
            lines = [line for line in f if line.strip() and not line.startswith("#")]
        input_contigs = {line.split("\t")[0] for line in lines}

        # Parse the file
        with open(path) as f:
            result = parse_metamaps_krona(path, f, ncbi_ranks)

        # Validate completeness
        assert len(result) > 0, "Should produce output"
        assert len(result) == len(input_contigs), "All contigs should be in output"

        # Check for duplicates
        output_contigs = [a.contig_name for a in result]
        assert len(output_contigs) == len(set(output_contigs)), "No duplicate contigs"
        assert set(output_contigs) == input_contigs, "All input contigs in output"

        # Verify mix of classified and unclassified reads
        classified = [r for r in result if r.clades != []]
        unclassified = [r for r in result if r.clades == []]
        assert len(classified) > 0, "Should have classified reads"
        assert len(unclassified) > 0, "Should have unclassified reads"

        # Test specific known contigs from fixture
        by_name = {r.contig_name: r for r in result}

        # ctg420_3x_l has taxID=0 (unclassified)
        assert "ctg420_3x_l" in by_name
        assert by_name["ctg420_3x_l"].clades == []

        # ctg631_2x_l has taxID=559292 (classified)
        assert "ctg631_2x_l" in by_name
        assert len(by_name["ctg631_2x_l"].clades) > 0

    def test_output_format(self, ncbi_ranks, fixtures_dir):
        """Output can be converted to valid string format."""
        path = fixtures_dir / "metamaps_valid.tsv"
        with open(path) as f:
            annotations = parse_metamaps_krona(path, f, ncbi_ranks)

        generic_annotations = [ncbi_ranks.generic_annotation(a) for a in annotations]

        for generic in generic_annotations:
            output_str = generic.to_string()
            assert "\t" in output_str
            parts = output_str.split("\t")
            assert len(parts) == 2
            assert ";" not in parts[0], "No semicolons in contig names"


class TestMetaMapsErrorHandling:
    """Test that invalid input produces helpful error messages."""

    def test_wrong_field_count(self, ncbi_ranks, tmp_path):
        """Wrong number of fields produces clear error."""
        test_file = tmp_path / "wrong_fields.tsv"
        test_file.write_text("contig1\t559292\n")  # Only 2 fields

        with open(test_file) as f:
            with pytest.raises(ValueError, match="expected 3.*got 2"):
                parse_metamaps_krona(test_file, f, ncbi_ranks)

    def test_invalid_taxid_format(self, ncbi_ranks, tmp_path):
        """Non-integer taxid produces clear error."""
        test_file = tmp_path / "invalid_taxid.tsv"
        test_file.write_text("contig1\tnotanumber\t1\n")

        with open(test_file) as f:
            with pytest.raises(ValueError, match="could not parse.*integer"):
                parse_metamaps_krona(test_file, f, ncbi_ranks)

    def test_unknown_taxid(self, ncbi_ranks, tmp_path):
        """Unknown taxid produces clear error with contig name."""
        test_file = tmp_path / "unknown_taxid.tsv"
        test_file.write_text("mystery_contig\t999999999\t1\n")

        with open(test_file) as f:
            with pytest.raises(ValueError, match="mystery_contig.*999999999.*not present"):
                parse_metamaps_krona(test_file, f, ncbi_ranks)

    def test_empty_file(self, ncbi_ranks, tmp_path):
        """Empty file is handled gracefully."""
        test_file = tmp_path / "empty.tsv"
        test_file.write_text("")

        with open(test_file) as f:
            result = parse_metamaps_krona(test_file, f, ncbi_ranks)

        assert result == [], "Empty input should produce empty output"


class TestMetaMapsFormatSpecifics:
    """Test MetaMaps-specific format quirks."""

    def test_comments_are_skipped(self, ncbi_ranks, tmp_path):
        """Lines starting with # are skipped."""
        test_file = tmp_path / "with_comments.tsv"
        content = (
            "# This is a comment\n"
            "contig1\t559292\t1\n"
            "# Another comment\n"
            "contig2\t0\t0\n"
        )
        test_file.write_text(content)

        with open(test_file) as f:
            result = parse_metamaps_krona(test_file, f, ncbi_ranks)

        assert len(result) == 2
        assert result[0].contig_name == "contig1"
        assert result[1].contig_name == "contig2"

    def test_taxid_zero_means_unclassified(self, ncbi_ranks, tmp_path):
        """taxID=0 produces empty clades (unclassified)."""
        test_file = tmp_path / "mixed.tsv"
        content = (
            "classified\t559292\t1\n"
            "unclassified\t0\t0\n"
        )
        test_file.write_text(content)

        with open(test_file) as f:
            result = parse_metamaps_krona(test_file, f, ncbi_ranks)

        assert len(result) == 2
        assert len(result[0].clades) > 0  # Classified
        assert result[1].clades == []  # Unclassified
