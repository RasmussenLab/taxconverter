from pathlib import Path
from taxconverter.common import GenericAnnotation, NCBIIdentifier


def parse_mmseqs_tsv(path: Path) -> list[GenericAnnotation]:
    result: list[GenericAnnotation] = []
    with open(path) as lines:
        for line_number, line in enumerate(lines, 1):
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"In MMseqs TSV file at {path}, on line {line_number}, "
                    f"expected 9 tab-separated fields, but got {len(fields)}. "
                    "Make sure to run mmseqs taxonomy with --tax-lineage"
                )

            result.append(
                GenericAnnotation(
                    fields[0], [NCBIIdentifier(i) for i in fields[8].split(";")]
                )
            )

    return result
