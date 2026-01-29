"""Tests for Metabuli parser - validates user-observable behavior."""

import pytest
from taxconverter.metabuli import parse_metabuli_files


class TestMetabuliValidInput:
    """Test that valid Metabuli input produces valid output."""

    def test_parse_valid_files(self, fixtures_dir):
        """Parse valid metabuli files and validate all output properties."""
        class_path = fixtures_dir / "metabuli_valid_classification.tsv"
        report_path = fixtures_dir / "metabuli_valid_report.tsv"

        # Read classification to count contigs
        with open(class_path) as f:
            lines = [line for line in f if line.strip() and not line.startswith("#")]
        input_contigs = {line.split("\t")[1] for line in lines}

        # Parse both files
        result = parse_metabuli_files(class_path, report_path)

        # Validate completeness
        assert len(result) > 0, "Should produce output"
        assert len(result) == len(input_contigs), "All contigs should be in output"

        # Check for duplicates
        output_contigs = [a.contig_name for a in result]
        assert len(output_contigs) == len(set(output_contigs)), "No duplicate contigs"
        assert set(output_contigs) == input_contigs, "All input contigs in output"

        # Verify mix of classified and unclassified
        classified = [r for r in result if r.clades != []]
        unclassified = [r for r in result if r.clades == []]
        assert len(classified) > 0, "Should have classified reads"
        assert len(unclassified) > 0, "Should have unclassified reads"

        # Test specific known contig from fixture
        by_name = {r.contig_name: r for r in result}

        # S16C209268 is unclassified (is_classified=0)
        assert "S16C209268" in by_name
        assert by_name["S16C209268"].clades == []

        # S16C9202 is classified
        assert "S16C9202" in by_name
        assert len(by_name["S16C9202"].clades) > 0

    def test_output_format(self, fixtures_dir):
        """Output can be converted to valid string format."""
        class_path = fixtures_dir / "metabuli_valid_classification.tsv"
        report_path = fixtures_dir / "metabuli_valid_report.tsv"

        result = parse_metabuli_files(class_path, report_path)

        for annotation in result:
            output_str = annotation.to_string()
            assert "\t" in output_str
            parts = output_str.split("\t")
            assert len(parts) == 2
            assert ";" not in parts[0], "No semicolons in contig names"


class TestMetabuliErrorHandling:
    """Test that invalid input produces helpful error messages."""

    def test_wrong_field_count_classified(self, tmp_path, fixtures_dir):
        """Classified reads must have 7 fields."""
        report_path = fixtures_dir / "metabuli_valid_report.tsv"
        class_file = tmp_path / "wrong_fields.tsv"
        class_file.write_text("1\tcontig1\t198454\t2997\t0.85619\n")  # Only 5 fields

        with pytest.raises(ValueError, match="expected 7.*got 5"):
            parse_metabuli_files(class_file, report_path)

    def test_wrong_field_count_unclassified(self, tmp_path, fixtures_dir):
        """Unclassified reads must have 6 fields."""
        report_path = fixtures_dir / "metabuli_valid_report.tsv"
        class_file = tmp_path / "wrong_fields_unclass.tsv"
        class_file.write_text("0\tcontig1\t0\n")  # Only 3 fields

        with pytest.raises(ValueError, match="expected 7.*got 3"):
            parse_metabuli_files(class_file, report_path)

    def test_invalid_is_classified_flag(self, tmp_path, fixtures_dir):
        """is_classified must be 0 or 1."""
        report_path = fixtures_dir / "metabuli_valid_report.tsv"
        class_file = tmp_path / "invalid_flag.tsv"
        class_file.write_text("2\tcontig1\t198454\t2997\t0.85619\tspecies\t198454:984\n")

        with pytest.raises(ValueError, match="expected first column.*'0' or '1'"):
            parse_metabuli_files(class_file, report_path)

    def test_taxid_not_in_report(self, tmp_path, fixtures_dir):
        """TaxID in classification must exist in report."""
        report_path = fixtures_dir / "metabuli_valid_report.tsv"
        class_file = tmp_path / "unknown_taxid.tsv"
        class_file.write_text("1\tcontig1\t999999999\t2997\t0.85619\tspecies\t999999999:984\n")

        with pytest.raises(ValueError, match="999999999.*not present.*report"):
            parse_metabuli_files(class_file, report_path)

    def test_empty_classification_file(self, tmp_path, fixtures_dir):
        """Empty classification file is handled gracefully."""
        report_path = fixtures_dir / "metabuli_valid_report.tsv"
        class_file = tmp_path / "empty.tsv"
        class_file.write_text("")

        result = parse_metabuli_files(class_file, report_path)
        assert result == []


class TestMetabuliFormatSpecifics:
    """Test Metabuli-specific format quirks."""

    def test_optional_header_handling(self, tmp_path, fixtures_dir):
        """Classification file may or may not have header."""
        report_path = fixtures_dir / "metabuli_valid_report.tsv"

        # Test with header
        class_with_header = tmp_path / "with_header.tsv"
        content = (
            "#is_classified\tname\ttaxID\tquery_length\tscore\trank\ttaxID:match_count\n"
            "1\tcontig1\t198454\t2997\t0.85619\tspecies\t198454:984\n"
        )
        class_with_header.write_text(content)
        result_with = parse_metabuli_files(class_with_header, report_path)

        # Test without header
        class_no_header = tmp_path / "no_header.tsv"
        class_no_header.write_text("1\tcontig1\t198454\t2997\t0.85619\tspecies\t198454:984\n")
        result_without = parse_metabuli_files(class_no_header, report_path)

        # Both should produce same result
        assert len(result_with) == 1
        assert len(result_without) == 1
        assert result_with[0].contig_name == result_without[0].contig_name

    def test_classified_and_unclassified(self, tmp_path, fixtures_dir):
        """Mixed classified and unclassified reads."""
        report_path = fixtures_dir / "metabuli_valid_report.tsv"
        class_file = tmp_path / "mixed.tsv"
        content = (
            "1\tclassified\t198454\t2997\t0.85619\tspecies\t198454:984\n"
            "0\tunclassified\t0\t29496\t0\tno rank\n"
        )
        class_file.write_text(content)

        result = parse_metabuli_files(class_file, report_path)

        assert len(result) == 2
        assert len(result[0].clades) > 0  # Classified
        assert result[1].clades == []  # Unclassified
        assert result[1].to_string() == "unclassified\t"
