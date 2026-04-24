def get_mem(snakemake, out_unit="MiB", mem_overhead_factor=0):
    """
    Obtain requested memory (from resources) and return in given units.
    If no memory resources found, return a value equivalent to 205 MiB.
    """
    import math

    # Store memory in MiB

    if mem_mb := snakemake.resources.get("mem_gb", None):
        mem_mb *= 1024
    else:
        mem_mb = snakemake.resources.get("mem_mb", 205)

    # Apply memory overhead
    if not (0 <= mem_overhead_factor < 1):
        raise ValueError(
            f"mem_overhead_factor must be >= 0 and < 1, got {mem_overhead_factor}"
        )
    mem_mb = math.floor(mem_mb * (1 - mem_overhead_factor))
    # Return memory
    if out_unit == "B":
        return mem_mb * 1024 * 1024
    elif out_unit == "KiB":
        return mem_mb * 1024
    elif out_unit == "MiB":
        return mem_mb
    elif out_unit == "GiB":
        return mem_mb / 1024
    else:
        raise ValueError("invalid output unit. Only B, KiB, MiB and GiB supported.")


def list_arg(cmd):
    """Turn command into list."""
    return list(filter(None, cmd.replace("=", " ").split(" ")))


def get_arg(arg, cmd):
    """Return position of argument on command."""
    try:
        return list_arg(cmd).index(arg)
    except ValueError:
        return None


def is_arg(arg, cmd):
    """Check presence of argument on command."""
    return get_arg(arg, cmd) is not None


def get_format(path, ignore_compression=True):
    """Get file format from extension, ignoring common compressions on user request"""
    from pathlib import Path

    if not path:
        raise ValueError("Path cannot be empty")
    exts = [s.lower() for s in Path(path).suffixes]

    if not exts:
        raise ValueError("Path must have an extension.")

    compression_fmt = {
        # https://en.wikipedia.org/wiki/Gzip
        ".gz": "gzip",
        ".gzip": "gzip",
        ".tgz": "gzip",
        # https://en.wikipedia.org/wiki/BGZF
        ".bgz": "bgzip",
        ".bgzip": "bgzip",
        # https://en.wikipedia.org/wiki/Bzip2
        ".bz2": "bzip2",
        # https://en.wikipedia.org/wiki/XZ_Utils
        ".xz": "lzma",
        ".lzma": "lzma",
        # https://en.wikipedia.org/wiki/List_of_file_formats
        ".mgz": "mgzip",
        ".zz": "zlib",
        ".z": "zlib",
        # https://github.com/sstadick/crabz/blob/91e58e3bdaaaf9838c14b5734947d82f2453be26/src/main.rs#L24
        ".snappy": "snap",
        ".snap": "snap",
        ".sz": "snap",
    }
    bioinfo_fmt = {
        # https://en.wikipedia.org/wiki/FASTA_format
        ".fa": "fasta",
        ".fas": "fasta",
        ".fna": "fasta",
        ".ffn": "fasta",
        ".faa": "fasta",
        ".fasta": "fasta",
        ".mpfa": "fasta",
        ".frn": "fasta",
        # https://en.wikipedia.org/wiki/FASTQ_format#File_extension
        ".fq": "fastq",
        ".fastq": "fastq",
    }

    if ignore_compression and exts[-1] in compression_fmt.keys():
        if len(exts) < 2:
            raise ValueError(
                "Compressed path must include a base extension before "
                "the compression suffix, e.g., '.vcf.gz'."
            )
        ext = exts[-2]
    else:
        ext = compression_fmt.get(exts[-1], exts[-1])
    ext = bioinfo_fmt.get(ext, ext)
    return ext.lstrip(".")


def move_files(snakemake, mapping, cmd="mv -v"):
    """
    Build shell move commands for relocating tool-produced files to named outputs.

    mapping must be a dict of {out_tag: source_path}. The out_tag must resolve
    to a single file path in snakemake.output.

    Example:
        mapping = {"tsv": "/tmp/tmp98723489/results/out.tsv"}

        # In the wrapper, one shell per move operation:
        for file in move_files(snakemake, mapping):
            shell("{file} {log}")

        # In the wrapper, one shell command for all move statements:
        move_cmds = "; ".join(move_files(snakemake, mapping))
        shell("(main_wrapper_cmd [...]; {move_cmds}) {log}")
    """

    cmds = []
    for out_tag, tool_out_name in mapping.items():
        out_name = snakemake.output.get(out_tag, "")
        if not out_name:
            raise KeyError(
                f"The wrapper requires the named output: {out_tag}. Please provide this named output."
            )
        cmds.append(f"{cmd} '{tool_out_name}' '{out_name}'")

    return cmds
