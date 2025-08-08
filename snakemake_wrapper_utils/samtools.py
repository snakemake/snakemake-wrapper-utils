import sys
from os import path
from snakemake_wrapper_utils.snakemake import is_arg


def infer_out_format(file_name):
    out_name, out_ext = path.splitext(file_name)
    return out_ext[1:].upper()


def get_samtools_opts(
    snakemake,
    parse_threads=True,
    parse_ref=True,
    parse_regions=True,
    parse_write_index=True,
    parse_output=True,
    parse_output_format=True,
    param_name="extra",
):
    """Obtain samtools_opts from output, params, and handle resource definitions in resources."""
    samtools_opts = ""
    extra = snakemake.params.get(param_name, "")
    idx = snakemake.output.get("idx", "")

    ###############
    ### Threads ###
    ###############
    if parse_threads:
        if is_arg("-@", extra) or is_arg("--threads", extra):
            sys.exit(
                "You have specified number of threads (`-@/--threads`) in `params.extra`; please use `threads`."
            )

        if snakemake.threads > 1:
            samtools_opts += f" --threads {snakemake.threads - 1}"

    ######################
    ### Reference file ###
    ######################
    if parse_ref:
        if is_arg("--reference", extra):
            sys.exit(
                "You have specified reference file (`--reference`) in `params.extra`; this is automatically infered from `ref` input file."
            )

        if snakemake.input.get("ref"):
            samtools_opts += f" --reference {snakemake.input.ref}"

    ####################
    ### Regions file ###
    ####################
    if parse_regions:
        if is_arg("--region-file", extra) or is_arg("--regions-file", extra):
            sys.exit(
                "You have specified regions file (`--region[s]-file`) in `params.extra`; this is automatically infered from `regions` input file."
            )

        if snakemake.input.get("regions"):
            samtools_opts += f" --region-file {snakemake.input.regions}"

    ###################
    ### Write index ###
    ###################
    if parse_write_index:
        if is_arg("--write-index", extra):
            sys.exit(
                "You have specified writing index (`--write-index`) in `params.extra`; this is automatically infered from `idx` output file."
            )
        if idx:
            samtools_opts += " --write-index"

    ###################
    ### Output file ###
    ###################
    if parse_output:
        if is_arg("-o", extra):
            sys.exit(
                "You have specified output file (`-o`) in `params.extra`; this is automatically infered from the first output file."
            )
        samtools_opts += f" -o {snakemake.output[0]}"
        if idx:
            samtools_opts += f"##idx##{idx}"

    #####################
    ### Output format ###
    #####################
    if parse_output_format:
        if is_arg("-O", extra) or is_arg("--output-fmt", extra):
            sys.exit(
                "You have specified output format (`-O/--output-fmt`) in `params.extra`; this is automatically infered from output file extension."
            )
        out_ext = infer_out_format(snakemake.output[0])
        samtools_opts += f" --output-fmt {out_ext}"

    return samtools_opts
