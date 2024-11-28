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
        raise valueError("invalid output unit. Only KiB, MiB and GiB supported.")


def is_arg(arg, cmd):
    """Check command for the presence of argument."""
    return arg in cmd.replace("=", " ").split(" ")


def get_format(path):
    from pathlib import Path
    """Get file format from extension, ignoring common compressions."""
    if not path:
        raise ValueError("Path cannot be empty")
    exts = Path(path).suffixes
    if not exts:
        raise ValueError("Path must have an extension")
    if exts[-1] in (".gz", ".bgz", ".bz2"):
        ext = exts[-2]
    else:
        ext = exts[-1]

    if ext in (".fq", ".fastq"):
        return "fastq"
    elif ext in (".fa", ".fas", ".fna", ".fasta"):
        return "fasta"
    else:
        return ext.lstrip(".").lower()
