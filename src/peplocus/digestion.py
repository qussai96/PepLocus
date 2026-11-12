"""Fast cleavage-context checks shared by exact, splice and one-edit searches."""
from .sequence import reverse_complement, translate

AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
METHIONINE_REMOVAL = frozenset("ACGPSTV")
NON_PROLINE = AMINO_ACIDS - {"P"}
C_FOLLOWING = NON_PROLINE | {"*"}
CONTEXT_FIELDS = ["upstream_aa", "downstream_aa", "query_n_terminal_context", "query_c_terminal_context",
                  "reference_n_terminal_context", "reference_c_terminal_context", "query_tryptic", "reference_tryptic",
                  "n_cleavage_change", "c_cleavage_change", "context_evidence"]


def cleavage_context(peptide, upstream="", downstream=""):
    """Check both ends; M-based exceptions establish plausibility, not initiation."""
    previous = upstream[-1:]  # upstream holds up to two residues, in protein order
    following = downstream[:1]
    n = "nontryptic" if previous in AMINO_ACIDS | {"*"} else "unknown"
    if previous in {"K", "R"} and peptide[:1] in NON_PROLINE:
        n = "tryptic"
    elif peptide.startswith("M"):
        n = "possible_protein_start"
    elif previous == "M" and peptide[:1] in METHIONINE_REMOVAL:
        n = "possible_methionine_removal"
    c = "nontryptic" if following in AMINO_ACIDS | {"*"} else "unknown"
    if following == "*":
        c = "translated_stop"
    elif peptide[-1:] in {"K", "R"} and following in NON_PROLINE:
        c = "tryptic"
    return dict(context_evidence="genome", upstream_aa=upstream[-2:], downstream_aa=following,
                n_terminal_context=n, c_terminal_context=c,
                enzyme_compatible=n not in {"nontryptic","unknown"} and c not in {"nontryptic","unknown"})


def compare_contexts(query, reference, upstream="", downstream=""):
    q=cleavage_context(query,upstream,downstream)
    ref=cleavage_context(reference,upstream,downstream)
    def change(which):
        flank=upstream[-1:] if which=="n" else downstream[:1]
        if flank not in AMINO_ACIDS | {"*"}:return "unresolved"
        if flank=="*":return "no_peptide_bond"
        before=ref[which+"_terminal_context"]=="tryptic"
        after=q[which+"_terminal_context"]=="tryptic"
        return "unchanged" if before==after else "gained" if after else "lost"
    return dict(upstream_aa=upstream[-2:],downstream_aa=downstream[:1],
        query_n_terminal_context=q["n_terminal_context"],query_c_terminal_context=q["c_terminal_context"],
        reference_n_terminal_context=ref["n_terminal_context"],reference_c_terminal_context=ref["c_terminal_context"],
        query_tryptic=q["enzyme_compatible"],reference_tryptic=ref["enzyme_compatible"],
        n_cleavage_change=change("n"),c_cleavage_change=change("c"),context_evidence="genome")


def protein_context(protein, start, end, query=None):
    reference=protein[start:end]
    return compare_contexts(reference if query is None else query,reference,
                            protein[max(0,start-2):start],protein[end:end+1])


def genomic_context(reference, chrom, strand, blocks, peptide, table=1):
    """Only outer flanks are needed; internal split codons remain in the peptide."""
    left, right = blocks[0][0], blocks[-1][1]
    size = reference.get_reference_length(chrom)
    if strand == "+":
        upstream = reference.fetch(chrom,max(0,left-6),left)
        downstream = reference.fetch(chrom,right,min(size,right+3))
        # At a contig edge use complete codons nearest the peptide.
        upstream = upstream[len(upstream)%3:]
    else:
        upstream = reverse_complement(reference.fetch(chrom,right,min(size,right+6)))
        downstream = reverse_complement(reference.fetch(chrom,max(0,left-3),left))
        upstream = upstream[len(upstream)%3:]
    return compare_contexts(peptide,peptide,translate(upstream,table),translate(downstream,table))


def extend_path(path, lengths, upstream_nt=6, downstream_nt=3):
    """Add flanking codons only to outer path ends, preserving all introns."""
    chrom,strand,segments,support,motifs,source=path
    segments=list(segments)
    if strand=="+":
        a,b=segments[0];segments[0]=(max(0,a-upstream_nt),b)
        a,b=segments[-1];segments[-1]=(a,min(lengths[chrom],b+downstream_nt))
    else:
        a,b=segments[0];segments[0]=(a,min(lengths[chrom],b+upstream_nt))
        a,b=segments[-1];segments[-1]=(max(0,a-downstream_nt),b)
    return chrom,strand,tuple(segments),support,motifs,source
