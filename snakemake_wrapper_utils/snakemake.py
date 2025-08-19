def get_mem(snakemake, out_unit="MiB"):
    """
    Obtain requested memory (from resources) and return in given units.
    If no memory resources found, return a value equivalent to 205.
    """

    # Store memory in MiB

    if mem_mb := snakemake.resources.get("mem_gb", None):
        mem_mb *= 1024
    else:
        mem_mb = snakemake.resources.get("mem_mb", 205)

    if out_unit == "KiB":
        return mem_mb * 1024
    elif out_unit == "MiB":
        return mem_mb
    elif out_unit == "GiB":
        return mem_mb / 1024
    else:
        raise ValueError("invalid output unit. Only KiB, MiB and GiB supported.")


def is_arg(arg, cmd):
    """Check command for the presence of argument."""
    return arg in cmd.replace("=", " ").split(" ")


def get_format(path):
    from pathlib import Path

    """Get file format from extension, ignoring common compressions."""
    if not path:
        raise ValueError("Path cannot be empty")
    exts = [s.lower() for s in Path(path).suffixes]
    if not exts:
        raise ValueError("Path must have an extension")
    if exts[-1] in (".gz", ".bgz", ".bz2", ".xz"):
        if len(exts) < 2:
            raise ValueError(
                "Compressed path must include a base extension before the compression suffix, e.g., '.vcf.gz'."
            )
        ext = exts[-2]
    else:
        ext = exts[-1]

    ext = ext.lstrip(".")
    if ext in ("fq", "fastq"):
        return "fastq"
    elif ext in ("fa", "fas", "fna", "fasta"):
        return "fasta"
    else:
        return ext


def move_files(snakemake, mapping, cmd="mv -v"):
    """
    Move file(s) produced by the tool to the named Snakemake outputs.

    mapping must be a dict of {out_tag: source_path}.

    Example:
        mapping = {"tsv": "/tmp/tmp98723489/results/out.tsv"}

    This moves /tmp/tmp98723489/results/out.tsv to snakemake.output["tsv"],
    redirecting stdout/stderr to snakemake.log.
    """

    from snakemake.shell import shell

    log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

    for out_tag, tool_out_name in mapping.items():
        out_name = snakemake.output.get(out_tag, "")
        if out_name:
            shell("{cmd} {tool_out_name:q} {out_name:q} {log}")
        else:
            raise KeyError(
                f"The wrapper requires the named output: {out_tag}. Please provide this named output."
            )
