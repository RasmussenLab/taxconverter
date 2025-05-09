import csv
import os
import sys
import time
import argparse
import pandas as pd
from loguru import logger
import taxconverter
from pathlib import Path
from typing import Iterable
import itertools

parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parentdir)

NCBI_LINEAGE = os.path.join(parentdir, 'data/clades.tsv')

TAXA_LEVELS = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
CHILD_ID = 'child_id'
PARENT_ID = 'parent_id'
CHILD_RANK = 'child_rank'
NAME = 'name'
KEY_COL = 'key'
LINEAGE_COL = 'lineage'
SEQ_COL = 'sequences'

COMMAND_CENTRIFUGE = 'centrifuge'
COMMAND_KRAKEN = 'kraken2'
COMMAND_METABULI = 'metabuli'
COMMAND_METAMAPS = 'metamaps'
COMMAND_MMSEQS = 'mmseqs2'

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



def all_to_mmseqs(df: pd.DataFrame):
    df['nan'] = 0
    df[LINEAGE_COL] = df[LINEAGE_COL].fillna('')
    df['len'] = df[LINEAGE_COL].str.split(';').map(len)
    df['last'] = df[LINEAGE_COL].str.split(';').str[-1]
    df['rank'] = df['len'].map(lambda x: TAXA_LEVELS[x-1])
    df = df[[SEQ_COL, 'nan', 'rank', 'last', 'nan', 'nan', 'nan', 'nan', LINEAGE_COL]]
    df.columns = list(range(9))
    return df


def all_to_taxvamb(df: pd.DataFrame):
    df[LINEAGE_COL] = df[LINEAGE_COL].fillna('')
    df = df[[SEQ_COL, LINEAGE_COL]]
    df.columns = ['contigs', 'predictions']
    return df


def ncbi_lineage():
    begintime = time.time()
    logger.info("Loading NCBI lineage")
    df_ncbi = pd.read_csv(NCBI_LINEAGE, quoting=csv.QUOTE_NONE, sep='\t')
    map_child_parent = {k: v for k, v in zip(df_ncbi[CHILD_ID].astype(str), df_ncbi[PARENT_ID].astype(str))}
    elapsed = round(time.time() - begintime, 2)
    logger.info(f"Loaded NCBI lineage with {len(map_child_parent)} entries in {elapsed} seconds")
    return map_child_parent


def get_lineage(tax_id, map_child_parent):
    if tax_id == 0 or tax_id == 1:
        logger.info(f"No lineage for ID: {tax_id}")
        return ''
    if tax_id not in map_child_parent:
        logger.info(f"ID not found in NCBI lineage: {tax_id}")
        return ''
    lineage = []
    while tax_id in map_child_parent and tax_id != '1':
        lineage.append(tax_id)
        tax_id = map_child_parent[tax_id]
    if tax_id != '1':
        logger.info(f"Parent node warning: {tax_id} is the root, but does not stem from 1. Returning empty lineage.")
        return ''
    return ';'.join(lineage[::-1])

def add_format_arguments(subparser):
    subparser.add_argument('-m', '--mmseqs-format', dest="mmseqs", action='store_true', help="convert to MMSeqs2 format (if you are using Taxometer)")


def add_one_filepath_arguments(subparser):
    subparser.add_argument('-i', '--input', dest="input", metavar="", type=str, help="path to the taxonomy annotations")
    subparser.add_argument('-o', '--output', dest="output", metavar="", type=str, help="path to save the converted annotations")
    add_format_arguments(subparser)


def add_metabuli_arguments(subparser):
    subparser.add_argument('-c', '--input-clas', dest="clas", metavar="", type=str, help="path to the Metabuli classification file")
    subparser.add_argument('-r', '--input-report', dest="report", metavar="", type=str, help="path to the Metabuli report file")
    subparser.add_argument('-o', '--output', dest="output", metavar="", type=str, help="path to save the converted annotations")
    add_format_arguments(subparser)


def check_float(lineno: int, s: str) -> float:
    try:
        return float(s)
    except ValueError:
        err = ValueError(
            f'On line {lineno}, could not parse column 0 as float: "{s}", perhaps format is wrong'
        )
        raise err from None


def check_int(lineno: int, s: str) -> int:
    f = check_float(lineno, s)
    if not f.is_integer():
        raise ValueError(
            f"On line {lineno}, numerical value is not integer: {s}, perhaps format is wrong."
        )
    return int(f)


def clean_result(clade_to_lineage: dict[int, list[str]], start_time: float) -> dict[int, str]:
    # Check if all have the same root name - if so, we delete it, since it's not really canonical
    if len(clade_to_lineage) > 0:
        all_same = True
        root_name = next(iter(clade_to_lineage.values()))
        for v in clade_to_lineage.values():
            if v[0] != root_name:
                all_same = False
                break

        if all_same:
            for v in clade_to_lineage.values():
                v.pop(0)

    result = {k: ";".join(v) for (k, v) in clade_to_lineage.items()}
    elapsed = time.time() - start_time
    
    logger.info(
        f"Loaded Metabuli dataset-specific lineage with {len(clade_to_lineage)} "
        f"entries in {elapsed:.2f} seconds"

    )
    return result


def metabuli_lineage(metabuli_report: Path) -> dict[int, str]:
    # Avoid an indentation level
    with open(metabuli_report) as file:
        return metabuli_from_iter(file)

def metabuli_from_iter(lines: Iterable[str]) -> dict[int, str]:
    start_time = time.time()
    clade_to_lineage: dict[int, list[str]] = dict()

    # Cache of last seen ranks - we need this to build the full lineage, since only
    # the current node is listed on each line, we need to keep track of its descendants
    ranks: list[str] = [""] * len(METABULI_CANONICAL_RANKS)
    last_rank_index = -1  # placeholder

    # Metabuli commit a17debb (2025-05-08) added a header to the file. So, the format may
    # or may not have the header depending on the version of Metabuli used.
    # Furthermore, the order of fields in this file have changed in the past, so we must
    # be fairly strict with parsing this file to avoid creating nonsense
    header = next(iter(lines), None)
    if header is None:
        return clean_result({}, start_time)
    elif header.strip() == "#clade_proportion\tclade_count\ttaxon_count\trank\ttaxID\tname":
        line_delta = 2
        iterator = lines
    else:
        # Add the line we just obtained back to the iterator
        iterator = itertools.chain([header], lines)
        line_delta = 1

    for lineno_minus_delta, line in enumerate(iterator):
        line = line.rstrip()
        lineno = lineno_minus_delta + line_delta
        newlineno = lineno
        # If we see an empty line, we check the rest of the file has empty lines. If so, return,
        # if not, throw an error since the file is malformatted. This is to handle trailing newlines
        # which editors sometimes add.
        if not line:
            for line in lines:
                newlineno += 1
                if line.rstrip():
                    raise ValueError(
                        f"Found empty line on line {lineno}, then nonempty on line {newlineno}"
                    )
            return clean_result(clade_to_lineage, start_time)

        (clade_proportion, clade_conut, taxon_count, rank, tax_id_str, clade) = (
            line.split("\t")
        )
        rank_index = METABULI_CANONICAL_RANKS.get(rank)

        if rank_index is None:
            raise ValueError(f'Unknown rank: "{rank}"')

        # Each successive rank has two more leading spaces in clade name.
        if not (
            len(clade) > 2 * rank_index
            # If rank_index is zero, then there are no leading spaces, so the isspace check fails
            and (rank_index == 0 or clade[: 2 * rank_index].isspace())
            and not clade[2 * rank_index].isspace()
        ):
            raise ValueError(
                f"On line {lineno}, leading spaces in clade name does not match rank"
            )

        stripped_clade = clade[2 * rank_index :]
        del clade  # avoid accidentally referring to unstripped clades after this point

        if stripped_clade in clade_to_lineage:
            raise ValueError(f'Duplicate clade seen: "{stripped_clade}"')

        # May or may not appear in output, but we need to ignore this if it's there.
        if stripped_clade == "unclassified":
            continue

        if ";" in stripped_clade:
            raise ValueError(
                f'Semicolon cannot appear in clade name "{stripped_clade}"'
            )

        # Parse numeric fields - we do this for safety, to make it more likely that if the
        # order of the (unlabelled) columns switch, as they seem to have done in
        # earlier versions of Metabuli, an error is thrown
        check_float(lineno, clade_proportion)
        check_int(lineno, clade_conut)
        check_int(lineno, taxon_count)
        clade_id = check_int(lineno, tax_id_str)

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
            ranks[rank_index] = stripped_clade
            clade_to_lineage[clade_id] = ranks[: rank_index + 1]
        else:
            raise ValueError(
                f"On line {lineno}, clade {stripped_clade} skips one or more ranks, or the rows are out of order"
            )

        last_rank_index = rank_index

    return clean_result(clade_to_lineage, start_time)

def main():

    def convert_to_unified(is_mmseqs):
        def decorator(func):
            def wrapper(*arguments, **kwargs):
                df = func(*arguments, **kwargs)
                if is_mmseqs:
                    df_result = all_to_mmseqs(df)
                    df_result.to_csv(args.output, header=None, sep='\t', index=None)
                else:
                    df_result = all_to_taxvamb(df)
                    df_result.to_csv(args.output, sep='\t', index=None)
            return wrapper
        return decorator

    doc = f"""
    Version: {'.'.join([str(i) for i in taxconverter.__version__])}

    Convert outputs of Metabuli, Centrifuge and Kraken2 to the unified format. The format is "contigs\tpredictions" and is accepted by TaxVAMB tool.
    Important: for the older release of Taxometer that uses MMSeqs2-like files, use the --mmseqs-format flag.
    As a result, an explicit full lineage is avaliable with each sequence id, using GTDB identifiers for Metabuli and MMSeqs2, NCBI identifiers for Centrifuge and Kraken2."""
    parser = argparse.ArgumentParser(
        prog="taxconverter",
        description=doc,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
    )
    helpos = parser.add_argument_group(title="Help and version", description=None)
    helpos.add_argument("-h", "--help", help="print help and exit", action="help")
    helpos.add_argument(
        "--version",
        action="version",
        version=f'Taxconverter {".".join(map(str, taxconverter.__version__))}',
    )

    subparsers = parser.add_subparsers(dest="subcommand")
    metabuli = subparsers.add_parser(
        COMMAND_METABULI,
        help="""
        Metabuli format converter
        """,
    )
    add_metabuli_arguments(metabuli)

    kraken = subparsers.add_parser(
        COMMAND_KRAKEN,
        help="""
        Kraken2 format converter
        """,
    )
    add_one_filepath_arguments(kraken)

    centrifuge = subparsers.add_parser(
        COMMAND_CENTRIFUGE,
        help="""
        Centrifuge format converter
        """,
    )
    add_one_filepath_arguments(centrifuge)

    metamaps = subparsers.add_parser(
        COMMAND_METAMAPS,
        help="""
        MetaMaps format converter
        """,
    )
    add_one_filepath_arguments(metamaps)

    mmseqs = subparsers.add_parser(
        COMMAND_MMSEQS,
        help="""
        MMSeqs2 format converter
        """,
    )
    add_one_filepath_arguments(mmseqs)

    args = parser.parse_args()

    @convert_to_unified(args.mmseqs)
    def centrifuge_data(filepath: str):
        map_ncbi = ncbi_lineage()
        begintime = time.time()
        df_centrifuge = pd.read_csv(filepath, delimiter='\t', usecols=['readID', 'taxID'])
        df_centrifuge['taxID'] = df_centrifuge['taxID'].astype(str)
        df_centrifuge[LINEAGE_COL] = df_centrifuge['taxID'].map(lambda x: get_lineage(x, map_ncbi))
        df_centrifuge[SEQ_COL] = df_centrifuge['readID']
        elapsed = round(time.time() - begintime, 2)
        logger.info(f"Converted Centrifuge to TaxVAMB/Taxometer format in {elapsed} seconds")
        return df_centrifuge


    @convert_to_unified(args.mmseqs)
    def metabuli_data(
            filepath_clas: Path,
            filepath_report: Path,
    ):
        map_lineage = metabuli_lineage(filepath_report)

        begintime = time.time()

        # Check if there is a header in the classification file.
        with open(filepath_clas) as file:
            clas_header = next(file, None)

        if clas_header is not None:
            if clas_header.startswith("#is_classified\tname\ttaxID\t"):
                pd_class_header = 0
            else:
                pd_class_header = None
        else:
            return pd.DataFrame()

        df_clas = pd.read_csv(
            filepath_clas,
            delimiter="\t",
            header=pd_class_header,
            dtype={0: int, 1: str, 2: int, 3: int, 4: float, 5: str, 6: str},
            names = list(range(7)),
        )

        tax_ids = df_clas[2]

        # Verify all tax ids are present in lineage file - otherwise failures will result in
        # silently returning zero annotations
        missing_id = next(filter(lambda x: x not in map_lineage, tax_ids), None)
        if missing_id is not None:
            raise ValueError(
                f"Tax ID {missing_id} from classifier file missing from lineage file"
            )

        df_clas[LINEAGE_COL] = tax_ids.map(map_lineage)
        df_clas[LINEAGE_COL] = df_clas[LINEAGE_COL].replace("unclassified", "")
        df_clas[LINEAGE_COL] = df_clas[LINEAGE_COL].str.replace('root;', '', regex=False)
        df_clas[SEQ_COL] = df_clas[1]
        df_clas_tax = df_clas[[SEQ_COL, LINEAGE_COL]]
        elapsed = round(time.time() - begintime, 2)
        logger.info(f"Converted Metabuli to TaxVAMB/Taxometer format in {elapsed} seconds")
        return df_clas_tax


    @convert_to_unified(args.mmseqs)
    def metamaps_data(filepath: str):
        df_ncbi = ncbi_lineage()
        begintime = time.time()
        df_metamaps = pd.read_csv(filepath, header=None, usecols=[0,1], delimiter='\t')
        df_metamaps.columns = ['readID', 'taxID']
        df_metamaps['taxID'] = df_metamaps['taxID'].astype(str)
        df_metamaps = pd.merge(df_metamaps, df_ncbi, left_on='taxID', right_on='tax_id', how='left')
        df_metamaps[SEQ_COL] = df_metamaps['readID']
        elapsed = round(time.time() - begintime, 2)
        logger.info(f"Converted MetaMaps to TaxVAMB/Taxometer format in {elapsed} seconds")
        return df_metamaps


    @convert_to_unified(args.mmseqs)
    def mmseqs_data(filepath: str):
        begintime = time.time()
        df_mmseqs = pd.read_csv(filepath, header=None, delimiter='\t')
        df_mmseqs[SEQ_COL] = df_mmseqs[0]
        df_mmseqs[LINEAGE_COL] = df_mmseqs[8]
        elapsed = round(time.time() - begintime, 2)
        logger.info(f"Converted MMseqs2 to TaxVAMB/Taxometer format in {elapsed} seconds")
        return df_mmseqs
    

    @convert_to_unified(args.mmseqs)
    def kraken_data(filepath: str):
        map_ncbi = ncbi_lineage()
        begintime = time.time()
        df_kraken = pd.read_csv(filepath, header=None, usecols=[1,2], delimiter='\t')
        df_kraken.columns = ['readID', 'taxID']
        df_kraken['taxID'] = df_kraken['taxID'].astype(str)
        df_kraken[SEQ_COL] = df_kraken['readID']
        df_kraken[LINEAGE_COL] = df_kraken['taxID'].map(lambda x: get_lineage(x, map_ncbi))
        elapsed = round(time.time() - begintime, 2)
        logger.info(f"Converted Kraken2 to TaxVAMB/Taxometer format in {elapsed} seconds")
        return df_kraken

    if args.subcommand == COMMAND_CENTRIFUGE:
        centrifuge_data(args.input)
    elif args.subcommand == COMMAND_KRAKEN:
        kraken_data(args.input)
    elif args.subcommand == COMMAND_METABULI:
        metabuli_data(args.clas, args.report)
    elif args.subcommand == COMMAND_METAMAPS:
        metamaps_data(args.input)
    elif args.subcommand == COMMAND_MMSEQS:
        mmseqs_data(args.input)
    else:
        assert False
