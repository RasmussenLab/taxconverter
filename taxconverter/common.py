import dataclasses
from typing import Optional, NewType, Self, Iterator
from pathlib import Path
import gzip

CANONICAL_RANK_NAMES = {
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
}

NCBI_HEADER = "child_id\tchild_rank\tparent_id\tname"

# 0-6 are the indices of TAXA_LEVELS
Rank = NewType("Rank", int)

RANK_FROM_STRING = {
    "domain": Rank(0),
    "phylum": Rank(1),
    "class": Rank(2),
    "order": Rank(3),
    "family": Rank(4),
    "genus": Rank(5),
    "species": Rank(6),
}

NCBIID = NewType("NCBIID", int)


@dataclasses.dataclass(frozen=True, slots=True)
class NCBIIdentifier:
    # Must not contain a newline, tab or semicolon
    content: str

    @classmethod
    def try_from_field(cls, field: str) -> Optional[Self]:
        if ";" in field:
            return None
        return cls(field)


@dataclasses.dataclass(frozen=True, slots=True)
class UnvalidatedIntAnnotation:
    contig_name: str
    leaf_clade: Optional[int]


@dataclasses.dataclass(frozen=True, slots=True)
class NCBIAnnotation:
    contig_name: str
    # The clade names, in order.
    # By default, these are assumed to represent: TAXA_LEVELS
    # This may be truncated, e.g. of length 2 if only domain and phylum is known
    clades: list[NCBIID]

    @classmethod
    def new_unknown(cls: type[Self], contig_name: str) -> Self:
        return cls(contig_name, [])


class NCBIRanks:
    # {child_id, (parent_id, child_rank, child_name)}
    # child rank is None if the rank is not canonical, e.g. Infraorder.
    __slots__ = ["child_to_parent"]
    child_to_parent: dict[NCBIID, tuple[NCBIID, Optional[Rank], NCBIIdentifier]]

    def annotation_from_int(
        self, int_annotation: UnvalidatedIntAnnotation
    ) -> NCBIAnnotation:
        integer = int_annotation.leaf_clade
        if integer is None or integer == 1:
            # Universal ancestor has identifier 1, and no parent, so is not in the parent dict
            return NCBIAnnotation.new_unknown(int_annotation.contig_name)

        value = self.child_to_parent.get(NCBIID(integer), None)
        if value is None:
            raise ValueError(
                f"Contig {repr(int_annotation.contig_name)} was annotated with NCBI ID {integer}, "
                "but this ID is not present in input NCBI file"
            )

        (parent_id, child_rank, _) = value
        # The file keeps intermediate ranks like "infraorder", but the parent of such
        # an intermediate rank is always canonical
        if child_rank is None:
            child_id = parent_id
            (parent_id, child_rank, _) = self.child_to_parent[child_id]
            # This property was validated in the construction of NCBIRanks
            assert child_rank is not None
        else:
            child_id = NCBIID(integer)

        clades: list[NCBIID] = [child_id]

        # This is the number of steps to take until child is SUPERKINGDOM
        for _ in range(int(child_rank)):
            child_id = parent_id
            (parent_id, child_rank, _) = self.child_to_parent[child_id]
            clades.append(parent_id)

        clades.reverse()

        return NCBIAnnotation(int_annotation.contig_name, clades)

    @classmethod
    def from_file(cls: type[Self], path: Path) -> Self:
        if path.suffix == ".gz":
            with gzip.open(path, "rt") as file:
                return cls.from_ncbi_lines(f"at path {path}", file)
        else:
            with open(path, "rt") as file:
                return cls.from_ncbi_lines(f"at path {path}", file)

    # This file is expected to be the output of
    # https://github.com/RasmussenLab/misc_scripts/blob/master/parse_ncbi_tax.jl
    @classmethod
    def from_ncbi_lines(cls: type[Self], where: str, lines: Iterator[str]) -> Self:
        first_line = next(lines, None)
        if first_line is None:
            raise ValueError(f"NCBI file {where} is empty")

        first_line = first_line.rstrip()

        if first_line != "child_id\tchild_rank\tparent_id\tname":
            raise ValueError(
                f"In first line of NCBI file {where}, got wrong header.\n"
                f"Expected:\n{repr(NCBI_HEADER)}\nGot:{repr(first_line)}"
            )

        map: dict[NCBIID, tuple[NCBIID, Optional[Rank], NCBIIdentifier]] = {}

        # Count from 2 because we just read the header above
        for line_number, line in enumerate(lines, 2):
            fields = line.split("\t")
            if len(fields) != 4:
                raise ValueError(
                    f"In NCBI file {where}, on line {line_number}, expected 4 fields, got {len(fields)}"
                )

            (child_id_str, child_rank_str, parent_id_str, name) = fields
            ncbi_id = NCBIIdentifier.try_from_field(name)
            if ncbi_id is None:
                raise ValueError(
                    f"In NCBI file {where}, on line {line_number}, clade name contains semicolon"
                )

            try:
                child_id = NCBIID(int(child_id_str))
                parent_id = NCBIID(int(parent_id_str))
            except ValueError:
                err = ValueError(
                    f"In NCBI file {where}, on line {line_number}, parent or child ID cannot be parsed as integer"
                )
                raise err from None

            rank = RANK_FROM_STRING.get(child_rank_str, None)
            map[child_id] = (parent_id, rank, ncbi_id)

        # Validation:
        # 1. Every non-canonical child has a canonical parent
        for (child_id), (parent_id, rank, ncbi_id) in map.items():
            if rank is None and parent_id != 1:
                if map[parent_id][1] is None:
                    raise ValueError(
                        f"Invalid NCBI file: Child {child_id} has a parent that is not "
                        f"one of the canonical ranks: {', '.join(CANONICAL_RANK_NAMES)}"
                    )

        instance = super().__new__(cls)
        instance.child_to_parent = map
        return instance
