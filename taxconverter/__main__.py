import csv
import os
import sys
import time
import argparse
import pandas as pd
from loguru import logger
import taxconverter


parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parentdir)

NCBI_LINEAGE = os.path.join(parentdir, 'data/clades.tsv')
METABULI_LINEAGE = os.path.join(parentdir, 'data/metabuli_lineage.csv')

TAXA_LEVELS = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
CHILD_ID = 'child_id'
PARENT_ID = 'parent_id'
CHILD_RANK = 'child_rank'
NAME = 'name'
LINEAGE_COL = 'lineage'
SEQ_COL = 'sequences'

COMMAND_CENTRIFUGE = 'centrifuge'
COMMAND_KRAKEN = 'kraken2'
COMMAND_METABULI = 'metabuli'
COMMAND_METAMAPS = 'metamaps'
COMMAND_MMSEQS = 'mmseqs2'


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
    lineage = []
    while tax_id != '1':
        lineage.append(tax_id)
        tax_id = map_child_parent.get(tax_id, '1')
    return ';'.join(lineage[::-1])


def metabuli_lineage():
    begintime = time.time()
    logger.info("Loading Metabuli lineage")
    df = pd.read_csv(METABULI_LINEAGE)
    elapsed = round(time.time() - begintime, 2)
    logger.info(f"Loaded Metabuli lineage with {len(df)} entries in {elapsed} seconds")
    return df


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
            filepath_clas: str,
            filepath_report: str,
    ):
        df_lineage = metabuli_lineage()

        begintime = time.time()
        df_report = pd.read_csv(filepath_report, delimiter='\t', header=None)
        map_ids_metabuli = {k: v.strip() for k, v in zip(df_report[4], df_report[5])}
        
        df_clas = pd.read_csv(filepath_clas, delimiter='\t', header=None)
        df_clas['label'] = df_clas[2].map(map_ids_metabuli)

        df_clas_tax = pd.merge(df_clas, df_lineage, left_on='label', right_on='key', how='left')
        df_clas_tax[SEQ_COL] = df_clas_tax[1]
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
