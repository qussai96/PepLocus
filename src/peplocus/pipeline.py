"""Mapping orchestration; user inputs are always read-only."""
from collections import Counter, deque
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import json
import csv
import shutil
from importlib.metadata import version
import os
import logging
import multiprocessing
from pathlib import Path
import platform
import tempfile
import time

import pysam

from . import __version__
from .annotation import read_annotation
from .discovery import discover
from .io import prepare_reference, read_peptides
from .junctions import extract_bam, read_junctions, splice_paths, write_junctions
from .output import HitStore
from .search import chunk_tasks, init_worker, make_automaton, search_chunk, search_path

LOG = logging.getLogger("peplocus")


def bounded_results(executor, tasks, workers):
    iterator = iter(tasks)
    pending = deque()
    for _ in range(workers * 2):
        task = next(iterator, None)
        if task is not None:
            pending.append(executor.submit(search_chunk, task))
    while pending:
        yield pending.popleft().result()
        task = next(iterator, None)
        if task is not None:
            pending.append(executor.submit(search_chunk, task))


def map_peptides(args):
    """Prevent concurrent writers from corrupting one output directory."""
    from .cli import validate_arguments
    validate_arguments(args)
    for name in ("genome", "peptides", "gtf", "bam", "junctions"):
        value=getattr(args,name,None)
        if value and not Path(value).is_file():
            raise ValueError(f"input {name} is not a readable file: {value}")
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    lock = out / ".peplocus.lock"
    try:
        handle = lock.open("x")
    except FileExistsError as error:
        raise ValueError(f"output is locked: {lock}; choose another output or remove a stale lock only after confirming its process has stopped") from error
    try:
        with handle:
            handle.write(f"pid={os.getpid()}\n")
        return _map_peptides(args)
    finally:
        lock.unlink(missing_ok=True)


def _map_peptides(args):
    args.il_equivalent = True
    started = time.monotonic()
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    products = ["all_mappings.bed", "tryptic_mappings.bed", "mappings.igv_display.bed", "mappings.tsv",
                "peptides.tsv", "unmapped.tsv", "junctions.tsv", "rescue", "summary.json"]
    legacy_products = ["mappings.bed", "rescue_mappings.bed", "rescue_mappings.tsv", "rescue_summary.json"]
    if any((out / name).exists() for name in products+legacy_products):
        raise ValueError(f"output directory already contains PepLocus results: {out}; choose a new --out directory")
    timings = {}
    peptides, input_stats = read_peptides(args.peptides, args.peptide_column, args.min_length)
    LOG.info("Loaded %s unique peptides (%s records)", len(peptides), input_stats["input_records"])
    genome, index, lengths = prepare_reference(args.genome, out / ".cache")
    LOG.info("Assembly: %s sequences, %s bp", len(lengths), sum(lengths.values()))
    LOG.info("One search; all placements and a separately counted tryptic subset")
    max_length = max(map(len, peptides))
    automaton = make_automaton(peptides, args.il_equivalent)
    limited = set()
    stages = {}
    with tempfile.TemporaryDirectory(prefix=".work-", dir=out) as work:
        work = Path(work)
        store = HitStore(work / "hits.sqlite", peptides, args.il_equivalent)
        try:
            stage_start = time.monotonic()
            tasks = chunk_tasks(lengths, args.chunk_size, max_length)
            init_args = (genome, index, automaton, peptides, args.genetic_code, args.il_equivalent)
            done = 0
            report = time.monotonic()
            if args.threads == 1:
                init_worker(*init_args)
                try:
                    for hits in map(search_chunk, tasks):
                        store.add(hits)
                        done += 1
                        if time.monotonic() - report > 20:
                            LOG.info("Contiguous search: %s chunks completed", done)
                            report = time.monotonic()
                finally:
                    from . import search
                    search._WORKER[0].close()
                    search._WORKER = None
            else:
                context = multiprocessing.get_context()
                with ProcessPoolExecutor(args.threads, initializer=init_worker, initargs=init_args, mp_context=context) as executor:
                    for hits in bounded_results(executor, tasks, args.threads):
                        store.add(hits)
                        done += 1
                        if time.monotonic() - report > 20:
                            LOG.info("Contiguous search: %s chunks completed", done)
                            report = time.monotonic()
            timings["contiguous_seconds"] = round(time.monotonic()-stage_start, 3)
            stages["contiguous_mapped_peptides"] = len(store.mapped())
            LOG.info("Contiguous search mapped %s peptides in %.1f s", stages["contiguous_mapped_peptides"], timings["contiguous_seconds"])
            with pysam.FastaFile(genome, filepath_index=index) as fasta:
                junctions = []
                if args.bam:
                    stage_start = time.monotonic()
                    LOG.info("Extracting splice junctions from BAM (no BAM index required)")
                    junctions, bam_stats = extract_bam(args.bam, fasta, lengths, args.min_mapq, args.min_read_anchor, args.min_junction_support, args.min_intron, args.max_intron)
                    stages["bam"] = bam_stats
                    timings["bam_extraction_seconds"] = round(time.monotonic()-stage_start, 3)
                    LOG.info("Retained %s BAM junctions in %.1f s", len(junctions), timings["bam_extraction_seconds"])
                if args.gtf:
                    stage_start = time.monotonic()
                    annotation_junctions, annotation_stats = read_annotation(args.gtf, fasta, lengths)
                    junctions.extend(annotation_junctions)
                    stages["annotation"] = annotation_stats
                    timings["annotation_extraction_seconds"] = round(time.monotonic()-stage_start, 3)
                    LOG.info("Loaded %s annotation junctions", len(annotation_junctions))
                if args.junctions:
                    junctions.extend(read_junctions(args.junctions, fasta, lengths, args.min_junction_support))
                if args.bam or args.junctions or args.gtf:
                    write_junctions(work / "junctions.tsv", junctions)
                    stage_start = time.monotonic()
                    path_stats = Counter()
                    buffer = []
                    report = time.monotonic()
                    for path in splice_paths(junctions, lengths, 3*max_length+6, args.max_junctions, args.max_paths_per_junction, path_stats):
                        buffer.extend(search_path(fasta, path, automaton, peptides, args.genetic_code, args.il_equivalent))
                        if len(buffer) >= 5000:
                            store.add(buffer)
                            buffer.clear()
                        if time.monotonic()-report > 20:
                            LOG.info("Splice search: %s paths searched", path_stats["splice_paths_searched"])
                            report = time.monotonic()
                    store.add(buffer)
                    stages["splice_paths"] = dict(path_stats)
                    stages["mapped_after_supported_junctions"] = len(store.mapped())
                    timings["supported_splice_seconds"] = round(time.monotonic()-stage_start, 3)
                    LOG.info("Supported junction search: %s mapped peptides total", stages["mapped_after_supported_junctions"])
                if args.discover_junctions:
                    stage_start = time.monotonic()
                    mapped = store.mapped()
                    targets = range(len(peptides)) if args.discover_all else (pid for pid in range(len(peptides)) if pid not in mapped)
                    novel_stats = Counter()
                    buffer = []
                    LOG.info("Discovering single-intron canonical junctions using terminal %s-aa anchors", args.novel_anchor)
                    for hit in discover(fasta, lengths, peptides, targets, args.genetic_code, args.il_equivalent, args.novel_anchor, args.min_intron, args.max_intron, args.max_seed_hits, args.max_novel_pairs, args.chunk_size, novel_stats, limited, cache_mb=args.novel_cache_mb):
                        buffer.append(hit)
                        if len(buffer) >= 5000:
                            store.add(buffer)
                            buffer.clear()
                    store.add(buffer)
                    stages["novel_discovery"] = dict(novel_stats)
                    timings["novel_discovery_seconds"] = round(time.monotonic()-stage_start, 3)
                rescue_rows={}
                if args.rescue_one_edit:
                    from .rescue_pipeline import run_rescue
                    mapped_ids = store.mapped()
                    rescue_summary = run_rescue(
                        genome, index, lengths, peptides,
                        (pid for pid in range(len(peptides)) if pid not in mapped_ids),
                        work, junctions=junctions, min_length=args.rescue_min_length,
                        threads=args.threads, chunk_size=args.chunk_size, table=args.genetic_code,
                        il_equivalent=args.il_equivalent, max_junctions=args.max_junctions,
                        max_paths=args.max_paths_per_junction, verify=not args.no_verify,
                        igv_display=args.igv_display_bed)
                    stages["one_edit_rescue"] = rescue_summary
                    timings["one_edit_rescue_seconds"] = rescue_summary["elapsed_seconds"]
                    if not rescue_summary["results"]["rescued_peptides"]:
                        shutil.rmtree(work/"rescue")
                    with (work/"rescue_peptides.tsv").open() as handle:
                        rescue_rows={row["peptide"]:row for row in csv.DictReader(handle,delimiter="\t")}
                stage_start = time.monotonic()
                LOG.info("Writing results%s", " and verifying every genomic mapping" if not args.no_verify else "")
                output_stats = store.finalize(work, fasta, args.genetic_code, not args.no_verify, limited, igv_display=args.igv_display_bed, rescue_rows=rescue_rows)
                timings["output_and_verification_seconds"] = round(time.monotonic()-stage_start, 3)
            tryptic_output_stats=output_stats.pop("tryptic_results")
            summary = {
                "peplocus_version": __version__, "finished_utc": datetime.now(timezone.utc).isoformat(),
                "python_version": platform.python_version(),
                "dependencies": {name: version(name) for name in ("numpy", "biopython", "pysam", "pyahocorasick")},
                "parameters": {k:v for k,v in vars(args).items() if k != "func"},
                "inputs": {str(Path(p).resolve()): {"size": Path(p).stat().st_size, "mtime_ns": Path(p).stat().st_mtime_ns} for p in (args.genome,args.peptides,args.bam,args.junctions,args.gtf) if p},
                "bed_label_conventions": {
                    "all_mappings.bed": "all exact sequence-compatible placements; biological N-to-C labels",
                    "tryptic_mappings.bed": "exact subset compatible with query cleavage context; scores use subset counts",
                    **({"mappings.igv_display.bed": "display-only; negative-strand peptide labels reversed; not for sequence identity comparisons"} if args.igv_display_bed else {}),
                },
                "assembly": {"sequences": len(lengths), "bases": sum(lengths.values())},
                "peptides": input_stats, "stages": stages, "results": output_stats,"tryptic_results":tryptic_output_stats,
                "enzyme_policy": {
                    "mode": "dual_outputs",
                    "rescue_classification": "query peptide in reference flanks; reference cleavage differences annotated, never used to discard an alignment",
                    "N_terminal_exceptions": "peptide begins M, or upstream M and peptide begins A/C/G/P/S/T/V; candidate initiation only",
                    "C_terminal_exception": "immediately followed by an in-frame translated stop",
                    "missed_cleavages": "unlimited",
                    "unknown_flanks": "do not establish cleavage or terminal boundaries",
                    "specificity_scope": "all and tryptic placement multiplicities counted independently; any compatible supported context qualifies a locus"
                },
                "search_scope": {
                    "contiguous": "all six frames and contigs; exhaustive for selected genetic code and I/L policy; enzyme labels do not remove placements",
                    "supported_splicing": "bounded combinations of supplied/BAM/annotation junctions; edges need not be read-phased",
                    "max_junctions_per_supported_path": args.max_junctions,
                    "splice_path_cap_reached": bool(stages.get("splice_paths",{}).get("splice_path_roots_capped",0)),
                    "novel_splicing": "disabled" if not args.discover_junctions else "canonical single-intron, bounded terminal anchors and intron lengths; unmapped peptides only unless --discover-all",
                    "coordinate_verification": not args.no_verify,
                    "mapping_confidence": "sequence compatibility only; not a peptide-spectrum or biological-origin probability"
                },
                "timings": timings, "elapsed_seconds": round(time.monotonic()-started,3)
            }
            (work / "summary.json").write_text(json.dumps(summary, indent=2)+"\n")
            store.close()
            store = None
            # summary.json is the completion marker and is published last.
            for name in products:
                if (work / name).exists():
                    (work / name).replace(out / name)
        finally:
            if store is not None:
                store.close()
    rescued=stages.get("one_edit_rescue",{}).get("results",{}).get("rescued_peptides",0)
    LOG.info("Finished: %s/%s peptides have exact locations (%s tryptic); %s additional peptides have one-edit candidates; %.1f s",
             output_stats["mapped_peptides"],len(peptides),tryptic_output_stats["mapped_peptides"],rescued,summary["elapsed_seconds"])
    return summary
