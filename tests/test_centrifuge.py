"""Tests for Centrifuge parser - validates user-observable behavior."""

import pytest
from taxconverter.centrifuge import parse_centrifuge, CENTRIFUGE_HEADER


class TestCentrifugeValidInput:
    """Test that valid Centrifuge input produces valid output."""

    def test_parse_valid_file(self, ncbi_ranks, fixtures_dir):
        """Parse valid centrifuge file and validate all output properties."""
        path = fixtures_dir / "centrifuge_valid.tsv"
        with open(path) as f:
            lines = f.readlines()

        # Parse the file
        with open(path) as f:
            result = parse_centrifuge(path, f, ncbi_ranks)

        # Extract input contig names (skip header line)
        input_contigs = {line.split("\t")[0] for line in lines[1:] if line.strip()}

        # Validate completeness
        assert len(result) > 0, "Should produce output"
        assert len(result) == len(input_contigs), "All contigs should be in output"

        # Check for duplicates
        output_contigs = [a.contig_name for a in result]
        assert len(output_contigs) == len(set(output_contigs)), "No duplicate contigs"

        # Verify all input contigs appear in output
        assert set(output_contigs) == input_contigs, "All input contigs in output"

    def test_output_format(self, ncbi_ranks, fixtures_dir):
        """Output can be converted to valid string format."""
        path = fixtures_dir / "centrifuge_valid.tsv"
        with open(path) as f:
            annotations = parse_centrifuge(path, f, ncbi_ranks)

        generic_annotations = [ncbi_ranks.generic_annotation(a) for a in annotations]

        for generic in generic_annotations:
            output_str = generic.to_string()
            assert "\t" in output_str
            parts = output_str.split("\t")
            assert len(parts) == 2
            assert ";" not in parts[0], "No semicolons in contig names"

    def test_unclassified_handling(self, ncbi_ranks, tmp_path):
        """taxID=0 produces empty clades (unclassified)."""
        test_file = tmp_path / "unclassified.tsv"
        content = (
            f"{CENTRIFUGE_HEADER}\n"
            "unclassified\tcid|0\t0\t0\t0\t1000\t1000\t0\n"
        )
        test_file.write_text(content)

        with open(test_file) as f:
            result = parse_centrifuge(test_file, f, ncbi_ranks)

        assert len(result) == 1
        assert result[0].clades == []
        assert ncbi_ranks.generic_annotation(result[0]).to_string() == "unclassified\t"


class TestCentrifugeErrorHandling:
    """Test that invalid input produces helpful error messages."""

    def test_header_validation(self, ncbi_ranks, tmp_path):
        """Wrong or missing header produces clear error."""
        test_file = tmp_path / "wrong_header.tsv"
        test_file.write_text("wrong\theader\nS2C0\tcid|1510\t1510\t100\t0\t1000\t1000\t1\n")

        with open(test_file) as f:
            with pytest.raises(ValueError, match="expected header"):
                parse_centrifuge(test_file, f, ncbi_ranks)

    def test_wrong_field_count(self, ncbi_ranks, tmp_path):
        """Wrong number of fields produces clear error."""
        test_file = tmp_path / "wrong_fields.tsv"
        test_file.write_text(f"{CENTRIFUGE_HEADER}\nS2C0\tcid|1510\t1510\n")

        with open(test_file) as f:
            with pytest.raises(ValueError, match="expected 8.*got 3"):
                parse_centrifuge(test_file, f, ncbi_ranks)

    def test_invalid_taxid_format(self, ncbi_ranks, tmp_path):
        """Non-integer taxid produces clear error."""
        test_file = tmp_path / "invalid_taxid.tsv"
        test_file.write_text(
            f"{CENTRIFUGE_HEADER}\nS2C0\tcid|1510\tnotanumber\t100\t0\t1000\t1000\t1\n"
        )

        with open(test_file) as f:
            with pytest.raises(ValueError, match="could not parse taxid"):
                parse_centrifuge(test_file, f, ncbi_ranks)

    def test_unknown_taxid(self, ncbi_ranks, tmp_path):
        """Unknown taxid produces clear error with contig name."""
        test_file = tmp_path / "unknown_taxid.tsv"
        test_file.write_text(
            f"{CENTRIFUGE_HEADER}\n"
            "mystery_contig\tcid|999999999\t999999999\t100\t0\t1000\t1000\t1\n"
        )

        with open(test_file) as f:
            with pytest.raises(ValueError, match="mystery_contig.*999999999.*not present"):
                parse_centrifuge(test_file, f, ncbi_ranks)

    def test_empty_file_after_header(self, ncbi_ranks, tmp_path):
        """Empty file (only header) is handled gracefully."""
        test_file = tmp_path / "empty.tsv"
        test_file.write_text(f"{CENTRIFUGE_HEADER}\n")

        with open(test_file) as f:
            result = parse_centrifuge(test_file, f, ncbi_ranks)

        assert result == []
