from pathlib import Path
from typing import Iterable
from taxconverter.common import UnvalidatedIntAnnotation, NCBIRanks, NCBIAnnotation

# From Kraken documentation
# Each sequence classified by Kraken results in a single line of output. Output lines contain five tab-delimited fields; from left to right, they are:

#     "C"/"U": one letter code indicating that the sequence was either classified or unclassified.
#     The sequence ID, obtained from the FASTA/FASTQ header.
#     The taxonomy ID Kraken used to label the sequence; this is 0 if the sequence is unclassified.
#     The length of the sequence in bp.
#     A space-delimited list indicating the LCA mapping of each k-mer in the sequence. For example, "562:13 561:4 A:31 0:1 562:3" would indicate that:
#         the first 13 k-mers mapped to taxonomy ID #562
#         the next 4 k-mers mapped to taxonomy ID #561
#         the next 31 k-mers contained an ambiguous nucleotide
#         the next k-mer was not in the database
#         the last 3 k-mers mapped to taxonomy ID #562


def parse_kraken(
    path: Path, lines: Iterable[str], ncbi: NCBIRanks
) -> list[NCBIAnnotation]:
    annotations: list[NCBIAnnotation] = []
    for line_number, line in enumerate(lines, 1):
        int_annotation = parse_kraken_line(path, line_number, line)
        annotation = ncbi.annotation_from_int(int_annotation)
        annotations.append(annotation)

    return annotations


def parse_kraken_line(
    path: Path, line_number: int, line: str
) -> UnvalidatedIntAnnotation:
    fields = line.split("\t")
    if len(fields) != 5:
        raise ValueError(
            f"In Kraken file at {path} on line {line_number} expected 5 tab-delimited fields, got {len(fields)}"
        )

    (cu, sequence_id, taxonomy_id, _, _) = fields
    if cu == "U":
        return UnvalidatedIntAnnotation(sequence_id, None)
    elif cu != "C":
        raise ValueError(
            f"In Kraken file at {path} on line {line_number}, first column must be 'C' or 'U'"
        )

    try:
        annotation = int(taxonomy_id)
    except ValueError:
        err = ValueError(
            f"In Kraken file at {path} on line {line_number}, could not parse annotation as integer"
        )
        raise err from None

    return UnvalidatedIntAnnotation(sequence_id, annotation)
