from typing import Iterator, NewType
from pathlib import Path
from taxconverter.common import (
    UnvalidatedIntAnnotation,
    NCBIIdentifier,
    GenericAnnotation,
)
import itertools

METABULI_CANONICAL_RANKS = {
    "no rank": 0,
    "superkingdom": 1,
    "phylum": 2,
    "class": 3,
    "order": 4,
    "family": 5,
    "genus": 6,
    "species": 7,
    "subspecies": 8,
}

METABULI_HEADER = (
    "#is_classified\tname\ttaxID\tquery_length\tscore\trank\ttaxID:match_count"
)

MetabuliID = NewType("MetabuliID", int)


def parse_metabuli_files(
    classification_path: Path,
    report_path: Path,
) -> list[GenericAnnotation]:
    with open(report_path) as file:
        id_to_identifiers = parse_metabuli_report(report_path, file)

    with open(classification_path) as file:
        unvalidated = parse_metabuli_classifications(classification_path, file)

    result: list[GenericAnnotation] = []
    for uv in unvalidated:
        integer = uv.leaf_clade
        if integer is None:
            result.append(GenericAnnotation.new_unknown(uv.contig_name))
            continue

        mb_id = MetabuliID(integer)
        identifiers = id_to_identifiers.get(mb_id, None)
        if identifiers is None:
            raise ValueError(
                f"Metabuli classification file at {classification_path} "
                f"contains taxon id {integer}, which is not present "
                f"in Metabuli report at {report_path}. "
                "Therefore, the lineage cannot be determined."
            )

        result.append(GenericAnnotation(uv.contig_name, identifiers))

    return result


def parse_metabuli_classifications(
    path: Path,
    lines: Iterator[str],
) -> list[UnvalidatedIntAnnotation]:
    result: list[UnvalidatedIntAnnotation] = []
    err_prefix = f"In file {path}, "

    header = next(lines, None)
    if header is None:
        return []
    if header.rstrip() != METABULI_HEADER:
        rest_lines = itertools.chain([header], lines)
        line_start = 1
    else:
        rest_lines = lines
        line_start = 2

    # Start from 2 since we skipped the header
    for line_number, line in enumerate(rest_lines, line_start):
        fields = line.split("\t")
        if len(fields) != 7:
            raise ValueError(
                err_prefix + f"on line {line_number}, expected 7 tab-separated fields, "
                f"but got {len(fields)}."
            )

        is_classified = fields[0]
        identifier = fields[1]
        tax_id_str = fields[2]

        if is_classified == "0" or tax_id_str == "0":
            result.append(UnvalidatedIntAnnotation(identifier, None))
            continue
        elif is_classified != "1":
            raise ValueError(
                err_prefix + f"on line {line_number}, "
                "expected first column to contain '0' or '1'"
            )

        try:
            taxid = int(tax_id_str)
        except ValueError:
            err = ValueError(
                err_prefix + f"on line {line_number}, could not parse taxid as integer"
            )
            raise err from None

        result.append(UnvalidatedIntAnnotation(identifier, taxid))

    return result


# The purpose of parse_metabuli report is to avoid having to parse the
# full list of all NCBI clades.
# Since the report includes all ancestors of all clades in the classification,
# This is a much smaller file to read, and so much more efficient.


# It looks like this: (without the leading # )
# (and with some numbers made up)
# 33.73   77571   77571   no rank 0        unclassified
# 100.0000        201563  0       no rank 1       root
# 100.0000        201563  221     superkingdom    14524     d__Bacteria
# 31.9642 64428   0       phylum  63357       p__Bacillota
# 31.9642 64428   25      class   63358         c__Bacilli
# 31.1248 62736   28      order   92156           o__Lactobacillales
# 29.4032 59266   12      family  92157             f__Streptococcaceae
# 26.7628 53944   5173    genus   92158               g__Streptococcus
# 4.5961  9264    8172    species 152593                s__Streptococcus agalactiae
#                            [ more lines ...]
# 0.0005  1       1       subspecies      153982                  RS_GCF_001592425.1
# 4.1957  8457    8312    species 139444                s__Streptococcus pyogenes
# 0.0174  35      35      subspecies      139697                  RS_GCF_000018125.1
def parse_metabuli_report(
    path: Path,
    lines: Iterator[str],
) -> dict[MetabuliID, list[NCBIIdentifier]]:
    clade_to_lineage: dict[MetabuliID, list[NCBIIdentifier]] = dict()

    # Cache of last seen ranks - we need this to build the full lineage, since only
    # the current node is listed on each line, we need to keep track of its descendants
    ranks: list[NCBIIdentifier] = [NCBIIdentifier.sentinel()] * len(
        METABULI_CANONICAL_RANKS
    )
    last_rank_index = -1  # placeholder

    # Metabuli commit a17debb (2025-05-08) added a header to the file. So, the format may
    # or may not have the header depending on the version of Metabuli used.
    # Furthermore, the order of fields in this file have changed in the past, so we must
    # be fairly strict with parsing this file to avoid creating nonsense
    header = next(lines, None)
    if header is None:
        return {}
    elif (
        header.strip()
        == "#clade_proportion\tclade_count\ttaxon_count\trank\ttaxID\tname"
    ):
        line_start = 2
        iterator = lines
    else:
        # Add the line we just obtained back to the iterator
        iterator = itertools.chain([header], lines)
        line_start = 1

    for line_number, line in enumerate(iterator, line_start):
        line = line.rstrip()
        # If we see an empty line, we check the rest of the file has empty lines. If so, return,
        # if not, throw an error since the file is malformatted. This is to handle trailing newlines
        # which editors sometimes add.
        if not line:
            seek_line_number = line_number
            for line in lines:
                seek_line_number += 1
                if line.rstrip():
                    raise ValueError(
                        f"In Metabuli report at {path}, found empty line on line {line_number}, "
                        "then nonempty on line {seek_line_number}"
                    )

            delete_universal_root(clade_to_lineage)
            return clade_to_lineage

        (clade_proportion, clade_count, taxon_count, rank, tax_id_str, clade) = (
            line.split("\t")
        )
        rank_index = METABULI_CANONICAL_RANKS.get(rank)

        if rank_index is None:
            raise ValueError(f'In Metabuli, found unknown rank: "{rank}"')

        # Each successive rank has two more leading spaces in clade name.
        if not (
            len(clade) > 2 * rank_index
            # If rank_index is zero, then there are no leading spaces, so the isspace check fails
            and (rank_index == 0 or clade[: 2 * rank_index].isspace())
            and not clade[2 * rank_index].isspace()
        ):
            raise ValueError(
                f"In Metabuli report at {path}, on line {line_number}, leading spaces in clade name does not match rank"
            )

        clade_identifier = NCBIIdentifier.try_from_field(clade[2 * rank_index :])
        if clade_identifier is None:
            raise ValueError(
                f"In Metabuli report at {path}, on line {line_number}, clade name contains semicolon"
            )
        del clade  # avoid accidentally referring to unstripped clades after this point

        # Parse numeric fields.
        # We only need the id, but we parse the other ones too.
        # We do this for safety, to make it more likely that if the
        # order of the (unlabelled) columns switch, as they seem to have done in
        # earlier versions of Metabuli, an error is thrown
        try:
            id = MetabuliID(int(tax_id_str))
            float(clade_proportion)
            int(clade_count)
            int(taxon_count)
        except ValueError:
            err = ValueError(
                f"In Metabuli report at {path}, on line {line_number}, "
                "expected columns 1, 2, 3, and 5 to contain a float, int, int and int, respectively."
            )
            raise err from None

        if id in clade_to_lineage:
            raise ValueError(f'Duplicate taxon ID seen: "{id}"')

        # Check that the rows are in correct order. This algorithm used by this parser
        # relies on the rows being well-ordered such that children of a clade directly
        # follows their parent or siblings. Any missing line will mess that up.
        # So, we add a check here.
        if (
            rank_index <= last_rank_index
            or rank_index == last_rank_index + 1
            # Special case if the non-canonical "root" is ever dropped from format
            or (rank_index == 1 and last_rank_index == -1)
        ):
            ranks[rank_index] = clade_identifier
            clade_to_lineage[id] = ranks[: rank_index + 1]
        else:
            raise ValueError(
                f"In Metabuli report at {path}, on line {line_number}, clade {clade_identifier} "
                "skips one or more ranks, or the rows are out of order"
            )

        last_rank_index = rank_index

    delete_universal_root(clade_to_lineage)
    return clade_to_lineage


# We have this function because Metabuli adds a "root" with taxid 1.
# This is not necessary, and it would be nicer to start at domain level,
# such that the levels corresponds to the canonical taxonomic ranks.
# So, here, if all clades descend from "root", we delete it
def delete_universal_root(
    clade_to_lineage: dict[MetabuliID, list[NCBIIdentifier]],
) -> None:
    if len(clade_to_lineage) == 0:
        return None

    root_name = next(iter(clade_to_lineage.values()))[0]

    for clade_names in clade_to_lineage.values():
        if clade_names[0] != root_name:
            return None

    for clade_names in clade_to_lineage.values():
        clade_names.pop(0)

    return None
