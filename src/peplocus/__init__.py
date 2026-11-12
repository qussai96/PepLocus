"""PepLocus: peptide coordinates directly from genome assemblies."""
__version__ = "0.5.6"


def map_peptides(assembly=None, peptides=None, out=None, **options):
    """Map peptide-file paths with CLI-equivalent options; return summary dict."""
    from .api import map_peptides as run
    return run(assembly, peptides, out, **options)


__all__ = ["__version__", "map_peptides"]
