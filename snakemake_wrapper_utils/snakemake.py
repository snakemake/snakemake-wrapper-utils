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


def _get_unnamed_arguments(parameter_list):
    """
    Get unnamed arguments of snakemake.input or snakemake.output.
    Find it as the first key in the parameter_list.

    If no unnamed arguments are found, return an empty list.

    Run as:
        _get_unnamed_arguments(snakemake.input)

    """
    keys_with_positions = parameter_list._names
    if len(keys_with_positions) == 0:
        # all arguments are unnamed
        return parameter_list

    first_key = next(iter(keys_with_positions.items()))
    n_unnamed_arguments = first_key[1][0]

    # as for input arguments, either is a string or a list of strings
    if n_unnamed_arguments == 1:
        return parameter_list[0]
    else:
        return [parameter_list[i] for i in range(n_unnamed_arguments)]
