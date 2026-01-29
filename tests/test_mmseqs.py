"""Tests for MMseqs2 parser - validates user-observable behavior."""

import pytest
from taxconverter.mmseqs import parse_mmseqs_tsv


class TestMMseqsValidInput:
    """Test that valid MMseqs2 input produces valid output."""

    def test_parse_valid_file(self, fixtures_dir):
        """Parse valid mmseqs file and validate all output properties."""
        path = fixtures_dir / "mmseqs_valid.tsv"

        # Read input to count contigs
        with open(path) as f:
            lines = f.readlines()
        input_contigs = {line.split("\t")[0] for line in lines if line.strip()}

        # Parse the file
        result = parse_mmseqs_tsv(path)

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

        # S14C105738 should be classified with GTDB lineage
        assert "S14C105738" in by_name
        assert len(by_name["S14C105738"].clades) > 1
        assert by_name["S14C105738"].clades[0].content == "d_Bacteria"

        # S14C139854 should be unclassified
        assert "S14C139854" in by_name
        assert by_name["S14C139854"].clades == []

    def test_output_format(self, fixtures_dir):
        """Output can be converted to valid string format."""
        path = fixtures_dir / "mmseqs_valid.tsv"
        annotations = parse_mmseqs_tsv(path)

        for annotation in annotations:
            output_str = annotation.to_string()

            # Basic format validation
            assert "\t" in output_str, "Output should be tab-separated"
            parts = output_str.split("\t")
            assert len(parts) == 2, "Output should be contig\\tlineage"

            contig_name, lineage = parts
            assert contig_name, "Contig name should not be empty"
            assert ";" not in contig_name, "Semicolons would break format"

            # Validate no forbidden characters
            assert "\n" not in output_str.rstrip("\n")
            assert "\r" not in output_str

            # Lineage validation
            if lineage:
                for part in lineage.split(";"):
                    assert "\t" not in part
                    assert "\n" not in part


class TestMMseqsErrorHandling:
    """Test that invalid input produces helpful error messages."""

    def test_wrong_field_count(self, tmp_path):
        """Wrong number of fields produces clear error with line number."""
        test_file = tmp_path / "wrong_fields.tsv"
        content = (
            "contig1\t60\tgenus\tC\t3\t3\t3\t1.000\td_Bacteria\n"
            "contig2\t60\tgenus\n"  # Only 3 fields on line 2
        )
        test_file.write_text(content)

        with pytest.raises(ValueError, match="line 2.*expected 9.*got 3"):
            parse_mmseqs_tsv(test_file)

    def test_empty_file(self, tmp_path):
        """Empty file is handled gracefully."""
        test_file = tmp_path / "empty.tsv"
        test_file.write_text("")

        result = parse_mmseqs_tsv(test_file)
        assert result == [], "Empty input should produce empty output"


class TestMMseqsFormatSpecifics:
    """Test MMseqs-specific format quirks."""

    def test_preserves_contig_order_and_gtdb_format(self, tmp_path):
        """Output preserves order and handles GTDB identifiers correctly."""
        test_file = tmp_path / "ordered.tsv"
        content = (
            "Z_contig\t60\tgenus\tC\t3\t3\t3\t1.000\td_Bacteria;p_Firmicutes\n"
            "A_contig\t0\tno rank\tunclassified\t1\t0\t0\t0.000\t\n"
            "M_contig\t48\tgenus\tN\t2\t2\t2\t1.000\td_Bacteria;p_Proteobacteria\n"
        )
        test_file.write_text(content)

        result = parse_mmseqs_tsv(test_file)

        # Order preserved
        assert result[0].contig_name == "Z_contig"
        assert result[1].contig_name == "A_contig"
        assert result[2].contig_name == "M_contig"

        # Classified have GTDB lineages
        assert [i.content for i in result[0].clades] == ["d_Bacteria", "p_Firmicutes"]
        assert [i.content for i in result[2].clades] == [
            "d_Bacteria",
            "p_Proteobacteria",
        ]

        # Unclassified have empty lineage
        assert result[1].clades == []
        assert result[1].to_string() == "A_contig\t"

    def test_get_mmseqs_rank(self, fixtures_dir):
        """MMseqs format can determine rank from lineage length."""
        path = fixtures_dir / "mmseqs_valid.tsv"
        annotations = parse_mmseqs_tsv(path)

        valid_ranks = {
            "domain",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
            "subspecies",
            "no rank",
        }

        for annotation in annotations:
            rank = annotation.get_mmseqs_rank()
            assert rank in valid_ranks, f"Invalid rank: {rank}"
