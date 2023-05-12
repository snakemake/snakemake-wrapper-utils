import sys
from snakemake_wrapper_utils.snakemake import get_mem, is_arg


def java_mem_xmx_error(params_key):
    return f"You have provided `-Xmx` in params.{params_key}. For Java memory specifications, please only use resources.mem_mb (for total memory reserved for the rule) and `params.java_mem_overhead_mb` (to specify any required non-heap overhead that needs to be set aside before determining the `-Xmx` value)."


def get_java_opts(snakemake, java_mem_overhead_factor=0.2):
    """Obtain java_opts from params, and handle resource definitions in resources."""

    java_opts = snakemake.params.get("java_opts", "")
    extra = snakemake.params.get("extra", "")
    assert 0.0 <= java_mem_overhead_factor <= 1.0

    # Getting memory.
    if "-Xmx" in java_opts:
        sys.exit(java_mem_xmx_error("java_opts"))
    if "-Xmx" in extra:
        sys.exit(java_mem_xmx_error("extra"))
    java_opts += " -Xmx{}M".format(
        round(get_mem(snakemake) * (1.0 - java_mem_overhead_factor))
    )

    # Getting java temp directory
    if "java_temp" in snakemake.output.keys():
        sys.exit(
            "You have specified `output.java_temp`; please use `resources.tmpdir`."
        )
    if is_arg("-Djava.io.tmpdir", java_opts):
        sys.exit(
            "You have specified `-Djava.io.tmpdir` in `params.java_opts`; please use `resources.tmpdir`."
        )
    if is_arg("-Djava.io.tmpdir", extra):
        sys.exit(
            "You have specified `-Djava.io.tmpdir` in `params.extra`; please use `resources.tmpdir`."
        )
    java_opts += f" -Djava.io.tmpdir={snakemake.resources.tmpdir}"

    return java_opts
