from pathlib import Path
from typing import Iterable
from taxconverter.common import NCBIRanks, NCBIAnnotation, UnvalidatedIntAnnotation

# Krona specification from https://github.com/marbl/Krona/wiki/Importing-NCBI-Taxonomy-IDs


def parse_metamaps_krona(
    path: Path, lines: Iterable[str], ncbi: NCBIRanks
) -> list[NCBIAnnotation]:
    result: list[NCBIAnnotation] = []
    for line_number, line in enumerate(lines, 1):
        # Comments are skipped in krona format
        if line.startswith("#"):
            continue

        fields = line.split("\t")
        if len(fields) != 3:
            raise ValueError(
                f"In Krona-formatted file at {path}, one line {line_number}"
                f"expected 3 tab-separated columns, but got {len(fields)}"
            )

        (identifier, taxid_str, _) = fields

        try:
            taxid = int(taxid_str)
        except ValueError:
            err = ValueError(
                f"In Krona file at {path} on line {line_number}, could not parse annotation as integer"
            )
            raise ValueError(err) from None

        # Annoyingly, this is undocumented in the krona format,
        # but Metamaps encodes a missing annotation as taxid zero.
        if taxid == 0:
            uninvalidated = UnvalidatedIntAnnotation(identifier, None)
        else:
            uninvalidated = UnvalidatedIntAnnotation(identifier, taxid)

        result.append(ncbi.annotation_from_int(uninvalidated))

    return result
