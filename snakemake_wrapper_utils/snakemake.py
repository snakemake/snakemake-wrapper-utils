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
        raise ValueError(f"mem_overhead_factor must be >= 0 and < 1, got {mem_overhead_factor}")
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
        raise ValueError("Path must have an extension.")
    if exts[-1] in (".gz", ".bgz", ".bz2", ".xz"):
        if len(exts) < 2:
            raise ValueError(
                "Compressed path must include a base extension before the compression suffix, e.g., '.vcf.gz'."
            )
        ext = exts[-2]
    else:
        ext = exts[-1]

    if ext in (".fq", ".fastq"):
        return "fastq"
    elif ext in (".fa", ".fas", ".fna", ".ffn", ".faa", ".fasta", ".mpfa", ".frn"):
        # https://en.wikipedia.org/wiki/FASTA_format
        return "fasta"
    else:
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
