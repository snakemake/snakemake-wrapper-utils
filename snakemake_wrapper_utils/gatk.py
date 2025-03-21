import sys
from snakemake_wrapper_utils.snakemake import is_arg


def get_gatk_opts(
    snakemake,
    parse_arg_file=False,
    parse_bam_index=False,
    parse_bam_md5=False,
    parse_vcf_index=False,
    parse_vcf_md5=False,
    parse_ref=False,
    parse_ref_dict=False,
    param_name="extra",
):
    """Obtain gatk_opts from input, output, params"""

    gatk_opts = ""
    extra = snakemake.params.get(param_name, "")

    ##########################
    ### Configuration file ###
    ##########################

    if parse_arg_file:
        if is_arg("--arguments_file", extra):
            sys.exit(
                "You have specified an argument file (`--argument_file`) in `params.extra`; this is automatically inferred from `input.arg_file`."
            )

        # Multiple argument files can be provided. Order matters.
        arg_file = snakemake.input.get("arg_file", "")
        if arg_file:
            if isinstance(arg_file, list):
                arg_file = " --argument_file ".join(arg_file)

            gatk_opts += f" --argument_file {arg_file}"

    ######################
    ### Reference file ###
    ######################

    if parse_ref:
        if is_arg("-R", extra) or is_arg("--reference", extra):
            sys.exit(
                "You have specified reference file (`-R,--reference`) in `params.extra`; this is automatically inferred from `input.ref`."
            )
        ref = snakemake.input.get("ref")
        if ref:
            gatk_opts += f" --reference {ref}"

    if parse_ref_dict:
        if is_arg("--sequence-dictionary", extra):
            sys.exit(
                "You have specified reference sequence dictionary (`--sequence-dictionary`) in `params.extra`; this is automatically inferred from `input.dict`."
            )
        dict = snakemake.input.get("dict", "")
        if dict:
            gatk_opts += f" --sequence-dictionary {dict}"

    ###########################
    ### Optional BAM output ###
    ###########################
    if parse_bam_index:
        if is_arg("--create-output-bam-index", extra) or is_arg("-OBI", extra):
            sys.exit(
                "You have specified bam index creation (`-OBI,--create-output-bam-index`) in `params.extra`; this is automatically inferred from `output.bam_bai`."
            )
        if snakemake.output.get("bam_bai"):
            gatk_opts += " --create-output-bam-index"

    if parse_bam_md5:
        if is_arg("--create-output-bam-md5", extra) or is_arg("-OBM", extra):
            sys.exit(
                "You have specified bam MD5-sum creation (`-OBM,--create-output-bam-md5`) in `params.extra`; this is automatically inferred from `output.bam_md5`."
            )
        if snakemake.output.get("bam_md5"):
            gatk_opts += " --create-output-bam-md5"

    ###########################
    ### Optional VCF output ###
    ###########################
    if parse_vcf_index:
        if is_arg("--create-output-variant-index", extra) or is_arg("-OVI", extra):
            sys.exit(
                "You have specified VCF index creation (`--OVI,--create-output-variant-index`) in `params.extra`; this is automatically inferred from `output.vcf_idx`."
            )
        if snakemake.output.get("vcf_idx"):
            gatk_opts += " --create-output-variant-index "

    if parse_vcf_md5:
        if is_arg("--create-output-variant-md5", extra) or is_arg("-OVM", extra):
            sys.exit(
                "You have specified VCF MD5-sum creation (`--OVI,--create-output-variant-index`) in `params.extra`; this is automatically inferred from `output.vcf_md5`."
            )
        if snakemake.output.get("vcf_md5"):
            gatk_opts += " --create-output-variant-md5 "

    return gatk_opts
