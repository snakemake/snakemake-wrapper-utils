def get_mem(snakemake, out_unit="MiB"):
    """
    Obtain requested memory (from resources) and return in given units.
    If no memory resources found, return a value equivalent to 205 .
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
    """Check command for the presence of argument"""
    if arg in cmd.replace("=", " ").split(" "):
        return True

    return False
