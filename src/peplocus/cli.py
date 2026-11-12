"""PepLocus command-line interface."""
import argparse
import csv
import json
import logging
from pathlib import Path
import sys
import tempfile

import pysam

from . import __version__
from .io import prepare_reference
from .junctions import extract_bam, write_junctions
from .annotation import read_annotation
from .output import EDIT_FIELDS
from .pipeline import map_peptides
from .rescue import one_edit
from .sequence import canonical, join_blocks, translate


class _HelpParser(argparse.ArgumentParser):
    """Keep tuning flags parseable while presenting a small everyday interface."""
    def _view(self, render):
        hidden=[] if getattr(self,"_show_advanced",False) else getattr(self,"_advanced_actions",[])
        saved=[(action,action.help) for action in hidden]
        try:
            for action,_ in saved:action.help=argparse.SUPPRESS
            return render()
        finally:
            for action,help_text in saved:action.help=help_text

    def parse_args(self, args=None, namespace=None):
        args=list(sys.argv[1:] if args is None else args)
        if args[:1]==["map"]:args=args[1:]  # compatibility with older scripts
        return super().parse_args(args, namespace)

    def format_help(self):return self._view(super().format_help)
    def format_usage(self):return self._view(super().format_usage)


class _AdvancedHelp(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        parser._show_advanced=True
        parser.print_help()
        parser.exit()


def positive(text):
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return value


def nonnegative(text):
    value = int(text)
    if value < 0:
        raise argparse.ArgumentTypeError("must be zero or a positive integer")
    return value


def junction_options(parser):
    parser.add_argument("--min-mapq", type=nonnegative, default=10, help="minimum BAM mapping quality (default: 10)")
    parser.add_argument("--min-read-anchor", type=positive, default=8, help="aligned bases on each side of a BAM intron (default: 8)")
    parser.add_argument("--min-junction-support", type=positive, default=2, help="minimum primary, nonduplicate alignment observations per junction (default: 2)")
    parser.add_argument("--min-intron", type=positive, default=20, help="minimum BAM/discovery intron length in bp (default: 20)")
    parser.add_argument("--max-intron", type=positive, default=500000, help="maximum BAM/discovery intron length in bp (default: 500000)")


def run_junctions(args):
    out = Path(args.out)
    if out.exists():
        raise ValueError(f"output already exists: {out}")
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="peplocus-junctions-") as cache:
        genome, index, lengths = prepare_reference(args.genome, cache)
        with pysam.FastaFile(genome, filepath_index=index) as fasta:
            junctions, stats = [], {}
            if args.bam:
                junctions, bam_stats = extract_bam(args.bam, fasta, lengths, args.min_mapq, args.min_read_anchor, args.min_junction_support, args.min_intron, args.max_intron)
                stats["bam"] = bam_stats
            if args.gtf:
                annotation_junctions, annotation_stats = read_annotation(args.gtf, fasta, lengths)
                junctions.extend(annotation_junctions)
                stats["annotation"] = annotation_stats
        write_junctions(out, junctions)
    print(json.dumps(stats, indent=2))


def run_verify(args):
    count = 0
    with tempfile.TemporaryDirectory(prefix="peplocus-verify-") as cache:
        genome, index, lengths = prepare_reference(args.genome, cache)
        with pysam.FastaFile(genome, filepath_index=index) as fasta, open(args.mappings) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not {"peptide", "chrom", "strand", "blocks", "matched_sequence"} <= set(reader.fieldnames or []):
                raise ValueError("expected a PepLocus mappings.tsv file")
            for row in reader:
                if not all(isinstance(row.get(name),str) and row[name] for name in ("peptide","chrom","strand","blocks","matched_sequence")):
                    raise ValueError("mapping row has empty required fields")
                blocks = json.loads(row["blocks"])
                if not isinstance(blocks,list) or not all(isinstance(block,list) and len(block)==2 and
                        all(type(value) is int for value in block) for block in blocks):
                    raise ValueError("mapping blocks must be a list of integer coordinate pairs")
                if row["strand"] not in {"+", "-"} or row["chrom"] not in lengths or not blocks:
                    raise ValueError(f"invalid mapping: {row}")
                previous = -1
                for a,b in blocks:
                    if not 0 <= a < b <= lengths[row["chrom"]] or a < previous:
                        raise ValueError(f"invalid genomic blocks: {blocks}")
                    previous = b
                if sum(b-a for a,b in blocks) != 3*len(row["matched_sequence"]):
                    raise ValueError("mapping blocks do not span complete matched codons")
                if "start" in row and (int(row["start"])!=blocks[0][0] or int(row["end"])!=blocks[-1][1]):
                    raise ValueError("mapping start/end disagree with blocks")
                observed = translate(join_blocks(fasta, row["chrom"], row["strand"], blocks), args.genetic_code)
                if row.get("match_type") == "one_edit":
                    edit = one_edit(row["peptide"], observed, args.il_equivalent)
                    valid = (observed == row["matched_sequence"] and edit is not None
                             and row.get("edit_distance") == "1"
                             and all(str(edit[key]) == row.get(key) for key in
                                     EDIT_FIELDS))
                else:
                    valid = observed == row["matched_sequence"] and canonical(observed,args.il_equivalent) == canonical(row["peptide"],args.il_equivalent)
                if not valid:
                    raise ValueError(f"mapping does not reproduce peptide/alignment: {row['peptide']} at {row['chrom']}:{blocks}")
                count += 1
    print(json.dumps({"verified_loci": count, "valid": True}))


def parser(command=None):
    if command in {"junctions", "verify"}:
        return utility_parser(command)
    mapping = _HelpParser(prog="peplocus",
        description=None,
        epilog="For tuning controls: peplocus --help-advanced",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    mapping.add_argument("--version", action="version", version=f"PepLocus {__version__}")
    mapping.add_argument("--help-advanced", action=_AdvancedHelp, nargs=0,
                         help="show advanced controls and exit")
    mapping.add_argument("-a", "--assembly", "--genome", dest="genome", metavar="ASSEMBLY", required=True, help="assembly FASTA (plain or gzip)")
    mapping.add_argument("-p", "--peptides", required=True, help="peptide list, FASTA, CSV or TSV")
    mapping.add_argument("--peptide-column", help="column name in CSV/TSV input")
    mapping.add_argument("--out", required=True, help="output directory without existing PepLocus results")
    mapping.add_argument("-g", "--gtf", "--gff", dest="gtf", help="annotation supplying splice junctions (GTF/GFF)")
    mapping.add_argument("-b", "--bam", help="optional RNA-seq BAM; no index required")
    mapping.add_argument("--junctions", help="optional reusable PepLocus junction TSV")
    mapping.add_argument("--threads", type=positive, default=8, help="worker processes (default: 8)")
    mapping.add_argument("--chunk-size", type=positive, default=2000000, help="genomic bases per worker task (default: 2000000)")
    mapping.add_argument("--min-length", type=positive, default=1, help="minimum peptide length; short matches can be highly ambiguous")
    mapping.add_argument("--il-equivalent", action="store_true", default=True, help=argparse.SUPPRESS)
    mapping.add_argument("--genetic-code", type=positive, default=1, help="NCBI translation table; applies to every contig (default: 1)")
    mapping.add_argument("--rescue-one-edit", action="store_true", help="rescue exact-unmapped peptides with one amino-acid edit")
    mapping.add_argument("--min-rescue-length", dest="rescue_min_length", metavar="MIN_LENGTH", type=positive, default=9, help="minimum query length for one-edit rescue (default: 9)")
    mapping.add_argument("--igv-display-bed", action="store_true", help="also write a display-only BED with reversed negative-strand labels and GAPMS-style colors")
    mapping.add_argument("--no-verify", action="store_true", help="skip final extraction/translation verification (enabled by default)")
    junction_options(mapping)
    mapping.add_argument("--max-junctions", type=positive, default=3, help="maximum introns crossed by a supported peptide path (default: 3)")
    mapping.add_argument("--max-paths-per-junction", type=nonnegative, default=1000, help="splice graph path cap per first intron; 0 means unlimited")
    mapping.add_argument("--discover-junctions", action="store_true", help="discover bounded single-intron alignments without annotation")
    mapping.add_argument("--discover-all", action="store_true", help="run discovery for mapped peptides too (requires --discover-junctions)")
    mapping.add_argument("--novel-anchor", type=positive, default=4, help="amino acids required at each terminal discovery anchor (default: 4)")
    mapping.add_argument("--novel-cache-mb", type=nonnegative, default=256, help="memory budget for caching small contigs during discovery; 0 disables it (default: 256 MiB)")
    mapping.add_argument("--max-seed-hits", type=nonnegative, default=2000, help="discovery seed-hit cap per contig/strand; 0 means unlimited")
    mapping.add_argument("--max-novel-pairs", type=nonnegative, default=10000, help="discovery candidate-junction cap per peptide/contig/strand; 0 means unlimited")
    common={"help","version","help_advanced","genome","peptides","out","gtf","bam","threads",
            "rescue_min_length","rescue_one_edit","discover_junctions"}
    mapping._advanced_actions=[action for action in mapping._actions if action.dest not in common]
    advanced=mapping.add_argument_group("Advanced controls")
    for action in mapping._advanced_actions:
        for group in mapping._action_groups:
            if group is not advanced and action in group._group_actions:group._group_actions.remove(action)
        advanced._group_actions.append(action)
    mapping.set_defaults(func=map_peptides, command="map")
    return mapping


def utility_parser(command):
    """Legacy utilities, kept for existing scripts; mapping needs no subcommand."""
    junction = argparse.ArgumentParser(prog="peplocus junctions")
    junction.add_argument("-a", "--assembly", "--genome", dest="genome", metavar="ASSEMBLY", required=True)
    junction.add_argument("-b", "--bam", help="RNA-seq BAM")
    junction.add_argument("-g", "--gtf", "--gff", dest="gtf", help="GTF or Helixer GFF3-style annotation")
    junction.add_argument("--out", required=True, help="new junction TSV path")
    junction_options(junction)
    junction.set_defaults(func=run_junctions, command="junctions")
    verify = argparse.ArgumentParser(prog="peplocus verify")
    verify.add_argument("-a", "--assembly", "--genome", dest="genome", metavar="ASSEMBLY", required=True)
    verify.add_argument("--mappings", required=True)
    verify.add_argument("--il-equivalent", action="store_true", default=True, help=argparse.SUPPRESS)
    verify.add_argument("--genetic-code", type=positive, default=1)
    verify.set_defaults(func=run_verify, command="verify")
    return junction if command == "junctions" else verify


def validate_arguments(args):
    if args.command == "junctions" and not (args.bam or args.gtf):
        raise ValueError("junctions requires --bam and/or --gtf")
    if hasattr(args,"min_intron") and args.min_intron > args.max_intron:
        raise ValueError("--min-intron cannot exceed --max-intron")
    if getattr(args,"discover_all",False) and not args.discover_junctions:
        raise ValueError("--discover-all requires --discover-junctions")
    from Bio.Data import CodonTable
    if getattr(args,"genetic_code",1) not in CodonTable.unambiguous_dna_by_id:
        raise ValueError("--genetic-code must be a valid NCBI DNA translation-table number")
    code=CodonTable.unambiguous_dna_by_id[getattr(args,"genetic_code",1)]
    if set(code.stop_codons) & set(code.forward_table):
        raise ValueError("context-dependent stop codons are unsupported for --genetic-code; choose a table with unambiguous stop semantics")
    if getattr(args,"rescue_min_length",2) < 2:
        raise ValueError("--min-rescue-length must be at least 2")


def main(argv=None):
    argv=list(sys.argv[1:] if argv is None else argv)
    command=argv.pop(0) if argv[:1] and argv[0] in {"junctions", "verify"} else None
    root = parser(command)
    args = root.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    try:
        validate_arguments(args)
    except ValueError as error:
        root.error(str(error))
    try:
        args.func(args)
    except (ValueError, OSError, KeyError) as error:
        logging.getLogger("peplocus").error("%s", error)
        return 2
    except KeyboardInterrupt:
        logging.getLogger("peplocus").error("Interrupted; no completed result summary was published")
        return 130
    return 0
