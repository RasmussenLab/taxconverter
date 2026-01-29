"""Tests for Kraken2 parser - validates user-observable behavior."""

import pytest
from taxconverter.kraken import parse_kraken


class TestKrakenValidInput:
    """Test that valid Kraken2 input produces valid output."""

    def test_parse_valid_file(self, ncbi_ranks, fixtures_dir):
        """Parse valid kraken file and validate all output properties."""
        path = fixtures_dir / "kraken_valid.tsv"
        with open(path) as f:
            lines = f.readlines()

        # Parse the file
        with open(path) as f:
            result = parse_kraken(path, f, ncbi_ranks)

        # Extract input contig names (kraken has no header)
        input_contigs = {line.split("\t")[1] for line in lines if line.strip()}

        # Validate completeness
        assert len(result) > 0, "Should produce output"
        assert len(result) == len(input_contigs), "All contigs should be in output"

        # Check for duplicates
        output_contigs = [a.contig_name for a in result]
        assert len(output_contigs) == len(set(output_contigs)), "No duplicate contigs"
        assert set(output_contigs) == input_contigs, "All input contigs in output"

        # Verify mix of classified and unclassified (based on fixture content)
        classified = [r for r in result if r.clades != []]
        unclassified = [r for r in result if r.clades == []]
        assert len(classified) > 0, "Should have classified reads"
        # Note: Only test if fixture has unclassified reads

    def test_output_format_and_classification(self, ncbi_ranks, tmp_path):
        """Output format is valid for both classified and unclassified reads."""
        test_file = tmp_path / "mixed.tsv"
        content = (
            "C\tclassified\t1510\t1000\t1510:500\n"
            "U\tunclassified\t0\t1000\t0:500\n"
        )
        test_file.write_text(content)

        with open(test_file) as f:
            result = parse_kraken(test_file, f, ncbi_ranks)

        # Validate classification
        assert len(result) == 2
        assert len(result[0].clades) > 0, "C flag should produce clades"
        assert result[1].clades == [], "U flag should produce empty clades"

        # Validate output format
        generic_annotations = [ncbi_ranks.generic_annotation(a) for a in result]
        for generic in generic_annotations:
            output_str = generic.to_string()
            assert "\t" in output_str
            parts = output_str.split("\t")
            assert len(parts) == 2
            assert ";" not in parts[0], "No semicolons in contig names"


class TestKrakenErrorHandling:
    """Test that invalid input produces helpful error messages."""

    def test_wrong_field_count(self, ncbi_ranks, tmp_path):
        """Wrong number of fields produces clear error."""
        test_file = tmp_path / "wrong_fields.tsv"
        test_file.write_text("C\tcontig1\t1510\n")  # Only 3 fields instead of 5

        with open(test_file) as f:
            with pytest.raises(ValueError, match="expected 5.*got 3"):
                parse_kraken(test_file, f, ncbi_ranks)

    def test_invalid_classification_flag(self, ncbi_ranks, tmp_path):
        """Invalid classification flag (not C or U) produces error."""
        test_file = tmp_path / "invalid_flag.tsv"
        test_file.write_text("X\tcontig1\t1510\t1000\t1510:500\n")

        with open(test_file) as f:
            with pytest.raises(ValueError, match="must be.*C.*or.*U"):
                parse_kraken(test_file, f, ncbi_ranks)

    def test_invalid_taxid_format(self, ncbi_ranks, tmp_path):
        """Non-integer taxid produces clear error."""
        test_file = tmp_path / "invalid_taxid.tsv"
        test_file.write_text("C\tcontig1\tnotanumber\t1000\t1510:500\n")

        with open(test_file) as f:
            with pytest.raises(ValueError, match="could not parse.*integer"):
                parse_kraken(test_file, f, ncbi_ranks)

    def test_unknown_taxid(self, ncbi_ranks, tmp_path):
        """Unknown taxid produces clear error with contig name."""
        test_file = tmp_path / "unknown_taxid.tsv"
        test_file.write_text("C\tmystery_contig\t999999999\t1000\t999999999:500\n")

        with open(test_file) as f:
            with pytest.raises(
                ValueError, match="mystery_contig.*999999999.*not present"
            ):
                parse_kraken(test_file, f, ncbi_ranks)

    def test_empty_file(self, ncbi_ranks, tmp_path):
        """Empty file is handled gracefully."""
        test_file = tmp_path / "empty.tsv"
        test_file.write_text("")

        with open(test_file) as f:
            result = parse_kraken(test_file, f, ncbi_ranks)

        assert result == [], "Empty input should produce empty output"


class TestKrakenFormatSpecifics:
    """Test Kraken-specific format quirks."""

    def test_no_header_and_c_u_flags(self, ncbi_ranks, tmp_path):
        """Kraken has no header and uses C/U classification flags."""
        test_file = tmp_path / "no_header.tsv"
        content = (
            "C\tcontig1\t1510\t1000\t1510:500\n"
            "U\tcontig2\t0\t1000\t0:500\n"
            "C\tcontig3\t520\t1000\t520:500\n"
        )
        test_file.write_text(content)

        with open(test_file) as f:
            result = parse_kraken(test_file, f, ncbi_ranks)

        # Should parse successfully without header
        assert len(result) == 3
        assert result[0].contig_name == "contig1"
        assert result[1].contig_name == "contig2"
        assert result[2].contig_name == "contig3"

        # Check classification
        assert len(result[0].clades) > 0
        assert result[1].clades == []
        assert len(result[2].clades) > 0
