from loguru import logger
import taxconverter
from taxconverter.common import (
    GenericAnnotation,
    NCBIAnnotation,
    NCBIRanks,
    NCBIIdentifier,
)
from pathlib import Path
import argparse
from typing import Optional, TextIO
from contextlib import nullcontext
import sys

# This file is stored in the manifest, and is required for this package
# to work correctly. It should automatically be downloaded when the package is.
NCBI_LINEAGE_PATH = Path(__file__).absolute().parent.parent / "data" / "clades.tsv.gz"

# List of supported subcommands
COMMAND_CENTRIFUGE = "centrifuge"
COMMAND_KRAKEN = "kraken2"
COMMAND_METABULI = "metabuli"
COMMAND_METAMAPS = "metamaps"
COMMAND_MMSEQS = "mmseqs2"

CLI_DOCUMENTATION = f"""
Version: {".".join([str(i) for i in taxconverter.__version__])}

Convert outputs of Metabuli, Centrifuge and Kraken2 to the unified format. The format is "contigs\tpredictions" and is accepted by TaxVAMB tool.
Important: for the older release of Taxometer that uses MMSeqs2-like files, use the --mmseqs-format flag.
As a result, an explicit full lineage is avaliable with each sequence id, using GTDB identifiers for Metabuli and MMSeqs2, NCBI identifiers for Centrifuge and Kraken2."""


def format_log(record) -> str:
    colors = {"WARNING": "red", "INFO": "green", "DEBUG": "blue", "ERROR": "red"}
    L = colors.get(record["level"].name, "blue")
    T = "red" if record["level"].name in ("WARNING", "ERROR") else "cyan"
    message = "<red>{message}</red>" if T == "red" else "{message}"
    time = f"<{T}>{{time:YYYY-MM-DD HH:mm:ss.SSS}}</{T}>"
    level = f"<b><{L}>{{level:<7}}</{L}></b>"
    return f"{time} | {level} | {message}\n{{exception}}"


# Replace default stderr logger with the one specified in format log above
logger.remove()
logger.add(sys.stderr, format=format_log)


def load_ncbi() -> taxconverter.common.NCBIRanks:
    # Allow users to load from either compressed or noncompressed (latter for speed)
    noncompressed = NCBI_LINEAGE_PATH.with_suffix("")
    if noncompressed.is_file():
        path = noncompressed
    elif NCBI_LINEAGE_PATH.is_file():
        path = NCBI_LINEAGE_PATH
    else:
        raise FileNotFoundError(
            f"Could not find Taxconverter lineage path at {str(noncompressed)}, "
            "with or without a '.gz' suffix.\n"
            "This file ought to come automatically with the Taxconverter installation,"
            "and is needed to extract full lineages from NCBI identifiers.\n"
            "Please manually download the file at github.com/RasmussenLab/taxconverter "
            "under 'Releases', and place it in the path as given above."
        )

    logger.info("Loading NCBI lineages.")

    ncbi = taxconverter.common.NCBIRanks.from_file(path)
    logger.info("\tDone loading NCBI lineages")
    return ncbi


def convert_ncbi_annotations(
    ncbi: NCBIRanks, annotations: list[NCBIAnnotation]
) -> list[GenericAnnotation]:
    result: list[GenericAnnotation] = []
    for annotation in annotations:
        result.append(ncbi.generic_annotation(annotation))
    return result


def write_output(
    destination: Optional[Path],
    annotations: list[GenericAnnotation],
    unassigned_clade_name: NCBIIdentifier,
    is_mmseqs: bool,
):
    dst_str = "stdout" if destination is None else destination
    logger.info(f"Writing output to {dst_str}")
    # Allow using a with-statement to safely close the file, but when
    # destination is None, do not close.
    if destination is None:
        context = nullcontext(sys.stdout)
    else:
        context = open(destination, "w")

    with context as output:
        if is_mmseqs:
            write_mmseqs_output(output, annotations, unassigned_clade_name)
        else:
            write_vamb_output(output, annotations, unassigned_clade_name)

    logger.info("\tDone writing output")


# MMseqs output is 9 columns:
# Name, zeros, rank, last, zeros, zeros, zeros, zeros, lineage
def write_mmseqs_output(
    output: TextIO,
    annotations: list[GenericAnnotation],
    unassigned_clade_name: NCBIIdentifier,
):
    unassigned_list = [unassigned_clade_name.content]
    for annotation in annotations:
        clades = annotation.clades
        clades = (
            [i.content for i in annotation.clades]
            if annotation.clades
            else unassigned_list
        )
        lineage = ";".join(clades)
        rank = annotation.get_mmseqs_rank()
        print(
            f"{annotation.contig_name}\t0\t{rank}\t{clades[-1]}\t0\t0\t0\t0\t{lineage}",
            file=output,
        )


def write_vamb_output(
    output: TextIO,
    annotations: list[GenericAnnotation],
    unassigned_clade_name: NCBIIdentifier,
):
    print("contigs\tpredictions", file=output)
    for annotation in annotations:
        # If annotation is empty, we use the unassigned clade name
        if not annotation.clades:
            print(f"{annotation.contig_name}\t{unassigned_clade_name.content}", file=output)
        else:
            print(annotation.to_string(), file=output)


def add_output(subparser: argparse.ArgumentParser):
    subparser.add_argument(
        "-o",
        "--output",
        dest="output",
        metavar="",
        type=Path,
        help="path to TSV output file [stdout]",
    )

    subparser.add_argument(
        "--unassigned",
        dest="unassigned",
        metavar="",
        type=str,
        default="unknown",
        help="Clade name for unassigned contigs ['unknown']",
    )

    subparser.add_argument(
        "-m",
        "--mmseqs-format",
        dest="output_mmseqs",
        action="store_true",
        help="output in MMSeqs2-like format",
    )


def main():
    parser = argparse.ArgumentParser(
        prog="taxconverter",
        description=CLI_DOCUMENTATION,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
    )
    helpos = parser.add_argument_group(title="Help and version", description=None)
    helpos.add_argument("-h", "--help", help="print help and exit", action="help")
    helpos.add_argument(
        "--version",
        action="version",
        version=f"Taxconverter {'.'.join(map(str, taxconverter.__version__))}",
    )

    subparsers = parser.add_subparsers(dest="subcommand")
    metabuli = subparsers.add_parser(
        COMMAND_METABULI,
        help="""
        Convert from Metabuli classification file
        """,
    )
    metabuli.add_argument(
        "-c",
        "--input-clas",
        dest="classifications",
        metavar="",
        required=True,
        type=Path,
        help="path to Metabuli *_classifications.tsv (required)",
    )
    metabuli.add_argument(
        "-r",
        "--input-report",
        dest="report",
        metavar="",
        required=True,
        type=Path,
        help="path to Metabuli *_report.tsv (required)",
    )
    add_output(metabuli)

    kraken = subparsers.add_parser(
        COMMAND_KRAKEN,
        help="""
        Convert from Kraken2 format
        """,
    )
    kraken.add_argument(
        "-i",
        "--input",
        dest="input",
        metavar="",
        required=True,
        type=Path,
        help="path to Kraken TSV file (required)",
    )
    add_output(kraken)

    centrifuge = subparsers.add_parser(
        COMMAND_CENTRIFUGE,
        help="""
        Convert from Centrifuge output
        """,
    )
    centrifuge.add_argument(
        "-i",
        "--input",
        dest="input",
        metavar="",
        required=True,
        type=Path,
        help="path to classification file (required)",
    )
    add_output(centrifuge)

    metamaps = subparsers.add_parser(
        COMMAND_METAMAPS,
        help="""
        Convert from Metamaps Krona format
        """,
    )
    metamaps.add_argument(
        "-i",
        "--input",
        dest="input",
        metavar="",
        required=True,
        type=Path,
        help="path to Metamaps Krona file (required)",
    )
    add_output(metamaps)

    mmseqs = subparsers.add_parser(
        COMMAND_MMSEQS,
        help="""
        Convert from MMseqs2 TSV with semicolon-sep lineage
        """,
    )
    mmseqs.add_argument(
        "-i",
        "--input",
        dest="input",
        metavar="",
        required=True,
        type=Path,
        help="path to MMseqs2 TSV with semicolon-sep lineage (required)",
    )
    add_output(mmseqs)

    args = parser.parse_args()

    for char in "\t;\r\n":
        if char in args.unassigned:
            raise ValueError(
                "Argument --unassigned must not contain semicolon, newlines, \\r or tab."
            )
    unassigned = NCBIIdentifier(args.unassigned)

    if args.output is not None:
        if args.output.exists():
            raise FileExistsError(args.output)

        if not args.output.parent.is_dir():
            raise ValueError(
                f"Parent of output {args.output} is not an existing directory"
            )

    if args.subcommand == COMMAND_CENTRIFUGE:
        if not args.input.is_file():
            raise FileNotFoundError(f"Centrifuge input file at {args.input}")

        ncbi = load_ncbi()

        logger.info("Loading Centrifuge input")
        logger.info(f"\tPath: {args.input}")
        with open(args.input) as file:
            annotations = taxconverter.centrifuge.parse_centrifuge(
                args.input, file, ncbi
            )
        logger.info("\tDone parsing Centrifuge input")
        generic = convert_ncbi_annotations(ncbi, annotations)
        write_output(args.output, generic, unassigned, args.output_mmseqs)

    elif args.subcommand == COMMAND_KRAKEN:
        if not args.input.is_file():
            raise FileNotFoundError(f"Kraken input file at {args.input}")

        ncbi = load_ncbi()

        logger.info("Loading Kraken input")
        logger.info(f"\tPath: {args.input}")
        with open(args.input) as file:
            annotations = taxconverter.kraken.parse_kraken(args.input, file, ncbi)
        logger.info("\tDone parsing Kraken input")
        generic = convert_ncbi_annotations(ncbi, annotations)
        write_output(args.output, generic, unassigned, args.output_mmseqs)

    elif args.subcommand == COMMAND_METABULI:
        if not args.classifications.is_file():
            raise FileNotFoundError(
                f"Metabuli classifications file at {args.classifications}"
            )

        if not args.report.is_file():
            raise FileNotFoundError(f"Metabuli report file at {args.report}")

        logger.info("Loading Metabuli input")
        logger.info(f"\tClassification paths: {args.classifications}")
        logger.info(f"\tReport paths: {args.report}")
        annotations = taxconverter.metabuli.parse_metabuli_files(
            args.classifications, args.report
        )
        logger.info("\tDone parsing Metabuli input")

        write_output(args.output, annotations, unassigned, args.output_mmseqs)
    elif args.subcommand == COMMAND_METAMAPS:
        if not args.input.is_file():
            raise FileNotFoundError(f"Metamaps Krona input file at {args.input}")

        ncbi = load_ncbi()

        logger.info("Loading Metamaps Krona input")
        logger.info(f"\tPath: {args.input}")
        with open(args.input) as file:
            annotations = taxconverter.metamaps.parse_metamaps_krona(
                args.input, file, ncbi
            )
        logger.info("\tDone parsing Metamaps input")

        generic = convert_ncbi_annotations(ncbi, annotations)
        write_output(args.output, generic, unassigned, args.output_mmseqs)

    elif args.subcommand == COMMAND_MMSEQS:
        if not args.input.is_file():
            raise FileNotFoundError(f"MMseqs taxonomy TSV file at {args.input}")

        logger.info("MMseqs taxonomy TSV file")
        logger.info(f"\tPath: {args.input}")
        annotations = taxconverter.mmseqs.parse_mmseqs_tsv(args.input)
        logger.info("\tDone parsing MMseqs input")
        write_output(args.output, annotations, unassigned, args.output_mmseqs)

    else:
        assert False
