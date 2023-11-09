import os,sys
import warnings
from .snakemake import get_mem, is_arg


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

    # Setting java temp directory

    # Check if it is not set multiple times

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
    # Check if env variable exists
    if "JDK_JAVA_OPTIONS" not in  os.environ:
        os.environ["JDK_JAVA_OPTIONS"] = ""
    else:
        # if tempdir is already defined in env variable, warn user and remove it
        if "-Djava.io.tmpdir" in os.environ["JDK_JAVA_OPTIONS"]:

            warnings.warn(
                " `-Djava.io.tmpdir` is already defined in `JDK_JAVA_OPTIONS`; Overwrite it with snakemake.resources.tmpdir."
            )

            # reges to remove the -Djava.io.tmpdir=.*\s
            import re
            os.environ["JDK_JAVA_OPTIONS"] = re.sub(
                "-Djava.io.tmpdir=.*\s", "", os.environ["JDK_JAVA_OPTIONS"]
            )
        
    # Append tempdir defined in resources
    os.environ["JDK_JAVA_OPTIONS"] += f" -Djava.io.tmpdir={snakemake.resources.tmpdir}"

    # return java_opts with memory
    return java_opts
