from pathlib import Path
from typing import Iterator
from taxconverter.common import (
    UnvalidatedIntAnnotation,
    NCBIRanks,
    NCBIAnnotation,
)

# From Centrifuge documentation
"""

The following example shows classification assignments for a read.  The assignment output has 8 columns.

    readID    seqID   taxID score	   2ndBestScore	   hitLength	queryLength	numMatches
    1_1	      gi|4    9646  4225	   0		       80	80		1

    The first column is the read ID from a raw sequencing read (e.g., 1_1 in the example).
    The second column is the sequence ID of the genomic sequence, where the read is classified (e.g., gi|4).
    The third column is the taxonomic ID of the genomic sequence in the second column (e.g., 9646).
    [...]
"""

CENTRIFUGE_HEADER = (
    "readID\tseqID\ttaxID\tscore\t2ndBestScore\thitLength\tqueryLength\tnumMatches"
)


def parse_centrifuge(
    path: Path, lines: Iterator[str], ncbi: NCBIRanks
) -> list[NCBIAnnotation]:
    result: list[NCBIAnnotation] = []
    header = next(lines, None)
    err_prefix = f"In Centrifuge classification file at {path}, "
    if header is None:
        raise ValueError(err_prefix + "got no header")
    if header.rstrip() != CENTRIFUGE_HEADER:
        raise ValueError(
            err_prefix
            + f"expected header:\n{repr(CENTRIFUGE_HEADER)}\nBut got a different header"
        )

    # Start from 2 since we skipped the header
    for line_number, line in enumerate(lines, 2):
        fields = line.split("\t")
        if len(fields) != 8:
            raise ValueError(
                err_prefix + f"on line {line_number}, expected 8 tab-separated fields, "
                f"but got {len(fields)}."
            )

        identifier = fields[0]
        tax_id_str = fields[2]

        try:
            taxid = int(tax_id_str)
        except ValueError:
            err = ValueError(
                err_prefix + f"on line {line_number}, could not parse taxid as integer"
            )
            raise err from None

        if taxid == 0:
            uninvalidated = UnvalidatedIntAnnotation(identifier, None)
        else:
            uninvalidated = UnvalidatedIntAnnotation(identifier, taxid)

        result.append(ncbi.annotation_from_int(uninvalidated))

    return result
