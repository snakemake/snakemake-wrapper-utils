import sys

def java_mem_xmx_error(unit, params_key):
    return f"You have specified resources.mem_{unit} and provided `-Xmx` in params.{params_key}. For Java memory specifications, please only use resources.mem_mb (for total memory reserved for the rule) and params.java_mem_overhead_mb (to specify any required non-heap overhead that needs to be set aside before determining the -Xmx value)."

def get_java_opts(snakemake, java_mem_overhead_factor=0.2):
    """Obtain java_opts from params, and handle resource definitions in resources."""

    java_opts = snakemake.params.get("java_opts", "")
    extra = snakemake.params.get("extra", "")
    assert 0.0 <= java_mem_overhead_factor <= 1.0

    # Getting memory in megabytes, if java opts is not filled with -Xmx parameter
    # By doing so, backward compatibility is preserved
    if "mem_mb" in snakemake.resources.keys():
        if "-Xmx" in java_opts:
            sys.exit(
                java_mem_xmx_error("mb", "java_opts")
            )
        if "-Xmx" in extra:
            sys.exit(
                java_mem_xmx_error("mb", "extra")
            )
        java_opts += " -Xmx{}M".format( round( snakemake.resources["mem_mb"] * (1.0 - java_mem_overhead_factor) ) )


    # Getting memory in gigabytes, for user convenience. Please prefer the use
    # of mem_mb over mem_gb as advised in documentation.
    elif "mem_gb" in snakemake.resources.keys():
        if "-Xmx" in java_opts:
            sys.exit(
                java_mem_xmx_error("gb", "java_opts")
            )
        if "-Xmx" in extra:
            sys.exit(
                java_mem_xmx_error("gb", "extra")
            )
        java_opts += " -Xmx{}G".format( round( snakemake.resources["mem_gb"] * (1.0 - java_mem_overhead_factor) ) )


    # Getting java temp directory from output files list, if -Djava.io.tmpdir
    # is not provided in java parameters. By doing so, backward compatibility is
    # not broken.
    if "java_temp" in snakemake.output.keys():
        if "-Djava.io.tmpdir" in java_opts:
            sys.exit(
                "You have specified output.java_temp and provided `-Djava.io.tmpdir` in params.java_opts. Please choose the one you intended and remove the other specification."
            )
        if "-Djava.io.tmpdir" in extra:
            sys.exit(
                "You have specified output.java_temp and provided `-Djava.io.tmpdir` in params.extra. Please choose the one you intended and remove the other specification."
            )
        java_opts += " -Djava.io.tmpdir={}".format(snakemake.output["java_temp"])

    return java_opts
