"""Python entry point using exactly the same defaults as the command line."""
import argparse


def map_peptides(assembly=None, peptides=None, out=None, *, threads=8, **options):
    """Map a peptide file; return the run summary dictionary.

    Keyword options use CLI names with underscores, e.g. gtf=..., bam=...,
    rescue_one_edit=True, min_rescue_length=9, chunk_size=2000000.
    Invalid options raise ValueError. Inputs are paths, never modified in place.
    """
    if "genome" in options:
        if assembly is not None:
            raise ValueError("use assembly only; genome is a legacy alias")
        assembly=options.pop("genome")
    if assembly is None or peptides is None or out is None:
        raise ValueError("assembly, peptides and out are required")
    from .cli import parser, validate_arguments
    root=parser()
    def error(message):raise ValueError(message)
    root.error=error
    actions={action.dest:action for action in root._actions if action.option_strings}
    if "min_rescue_length" in options:
        if "rescue_min_length" in options:
            raise ValueError("use only min_rescue_length (rescue_min_length is a legacy alias)")
        options["rescue_min_length"]=options.pop("min_rescue_length")
    if "il_equivalent" in options and options["il_equivalent"] is not True:
        raise ValueError("I/L equivalence is always enabled; il_equivalent must be True")
    values=dict(options,threads=threads)
    argv=["--assembly",str(assembly),"--peptides",str(peptides),"--out",str(out)]
    for key,value in values.items():
        if key not in actions or key in {"help","help_advanced","version","genome","peptides","out"}:
            raise ValueError("unknown or reserved mapping option: "+key)
        action=actions[key]
        flag=action.option_strings[0]
        if isinstance(action,argparse._StoreTrueAction):
            if not isinstance(value,bool):raise ValueError(key+" must be boolean")
            if value:argv.append(flag)
        elif value is not None:
            argv.extend([flag,str(value)])
    args=root.parse_args(argv)
    validate_arguments(args)
    return args.func(args)
