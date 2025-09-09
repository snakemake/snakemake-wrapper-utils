import sys
from snakemake_wrapper_utils.snakemake import get_mem, is_arg


def infer_out_format(output, uncompressed_bcf=False):
    if output.endswith(".vcf"):
        return "v"
    elif output.endswith(".vcf.gz"):
        return "z"
    elif output.endswith(".bcf"):
        if uncompressed_bcf:
            return "u"
        else:
            return "b"
    else:
        raise ValueError("invalid output file format ('.vcf', '.vcf.gz', '.bcf').")


def get_bcftools_opts(
    snakemake,
    parse_threads=True,
    parse_ref=True,
    parse_regions=True,
    parse_samples=True,
    parse_targets=True,
    parse_output=True,
    parse_output_format=True,
    parse_memory=True,
):
    """Obtain bcftools_opts from output, params, and handle resource definitions."""
    bcftools_opts = ""
    extra = snakemake.params.get("extra", "")

    ###############
    ### Threads ###
    ###############
    if parse_threads:
        if is_arg("--threads", extra):
            sys.exit(
                "You have specified number of threads (`--threads`) in `params.extra`; please use `threads`."
            )
        bcftools_opts += (
            ""
            if snakemake.threads <= 1
            else "--threads {}".format(snakemake.threads - 1)
        )

    ######################
    ### Reference file ###
    ######################
    if parse_ref:
        if is_arg("-f", extra) or is_arg("--fasta-ref", extra):
            sys.exit(
                "You have specified reference file (`-f/--fasta-ref`) in `params.extra`; this is automatically infered from the `ref` input file."
            )

        if snakemake.input.get("ref"):
            bcftools_opts += f" --fasta-ref {snakemake.input.ref}"

    ####################
    ### Regions file ###
    ####################
    if parse_regions:
        if is_arg("--regions-file", extra) or is_arg("-R", extra):
            sys.exit(
                "You have specified regions file (`-R/--regions-file`) in `params.extra`; this is automatically infered from the `regions` input file."
            )

        if snakemake.input.get("regions"):
            if not snakemake.input.get("index"):
                sys.exit(
                    "You have specified a `regions` file in `input:`; this implies the `--regions-file` option of bcftools and thus also requires an `index=` file specified in the `input:`, but none was found."
                )
            bcftools_opts += f" --regions-file {snakemake.input.regions}"

    ####################
    ### Samples file ###
    ####################
    if parse_samples:
        if is_arg("-S", extra) or is_arg("--samples-file", extra):
            sys.exit(
                "You have specified samples file (`-S/--samples-file`) in `params.extra`; this is automatically infered from the `samples` input file."
            )

        if snakemake.input.get("samples"):
            bcftools_opts += f" --samples-file {snakemake.input.samples}"

    ####################
    ### Targets file ###
    ####################
    if parse_targets:
        if is_arg("-T", extra) or is_arg("--targets-file", extra):
            sys.exit(
                "You have specified targets file (`-T/--targets-file`) in `params.extra`; this is automatically infered from the `targets` input file."
            )

        if snakemake.input.get("targets"):
            bcftools_opts += f" --targets-file {snakemake.input.targets}"

    ###################
    ### Output file ###
    ###################
    if parse_output:
        if is_arg("-o", extra) or is_arg("--output", extra):
            sys.exit(
                "You have specified output file (`-o/--output`) in `params.extra`; this is automatically infered from the first output file."
            )
        bcftools_opts += f" --output {snakemake.output[0]}"

    #####################
    ### Output format ###
    #####################
    if parse_output_format:
        if is_arg("-O", extra) or is_arg("--output-type", extra):
            sys.exit(
                "You have specified output format (`-O/--output-type`) in `params.extra`; this is automatically infered from output file extension."
            )

        out_format = infer_out_format(
            snakemake.output[0], snakemake.params.get("uncompressed_bcf", False)
        )
        bcftools_opts += f" --output-type {out_format}"

    ##############
    ### Memory ###
    ##############
    if parse_memory:
        if is_arg("-m", extra) or is_arg("--max-mem", extra):
            sys.exit(
                "You have provided `-m/--max-mem` in `params.extra`; please use `resources.mem_mb`."
            )
        bcftools_opts += " --max-mem {}M".format(get_mem(snakemake))

    ################
    ### Temp dir ###
    ################
    if (
        is_arg("-T", extra)
        or is_arg("--temp-dir", extra)
        or is_arg("--temp-prefix", extra)
    ):
        sys.exit(
            "You have provided `-T/--temp-dir/--temp-prefix` in `params.extra`; please use `resources.tmpdir`."
        )

    return bcftools_opts
