# Changelog

## [0.6.0](https://www.github.com/snakemake/snakemake-wrapper-utils/compare/v0.5.3...v0.6.0) (2023-06-21)


### Features

* allow to configure the param name parsed by get_samtools_opts ([#30](https://www.github.com/snakemake/snakemake-wrapper-utils/issues/30)) ([84aa151](https://www.github.com/snakemake/snakemake-wrapper-utils/commit/84aa1515bc2ea34b702490b80b40e275ddc1e4af))


### Bug Fixes

* fix bug when short arguments are contained in long arguments ([#28](https://www.github.com/snakemake/snakemake-wrapper-utils/issues/28)) ([6b9b6db](https://www.github.com/snakemake/snakemake-wrapper-utils/commit/6b9b6db33ed406a88f1f1a83c7e1d1f5cf77d18e))

### [0.5.3](https://www.github.com/snakemake/snakemake-wrapper-utils/compare/v0.5.2...v0.5.3) (2023-03-14)


### Bug Fixes

* check for `input: index=""`, when `input: regions=""` is specified, as the `--regions-file` option requires an index ([#26](https://www.github.com/snakemake/snakemake-wrapper-utils/issues/26)) ([7e11735](https://www.github.com/snakemake/snakemake-wrapper-utils/commit/7e117351211369e4f58753845f3fc19d5fad7606))

### [0.5.2](https://www.github.com/snakemake/snakemake-wrapper-utils/compare/v0.5.1...v0.5.2) (2023-01-09)


### Bug Fixes

* introduce default 20% deduction from resources.mem_mb and override via params.java_mem_overhead_mv ([#24](https://www.github.com/snakemake/snakemake-wrapper-utils/issues/24)) ([71403bd](https://www.github.com/snakemake/snakemake-wrapper-utils/commit/71403bd4a843cc66bee28c6a11a279654b2b1857))

### [0.5.1](https://www.github.com/snakemake/snakemake-wrapper-utils/compare/v0.5.0...v0.5.1) (2022-12-09)


### Bug Fixes

* typo when checking for regions file in params.extra ([#22](https://www.github.com/snakemake/snakemake-wrapper-utils/issues/22)) ([fa70493](https://www.github.com/snakemake/snakemake-wrapper-utils/commit/fa704938530accc957f48b267dde051cbd8d20fc))

## [0.5.0](https://www.github.com/snakemake/snakemake-wrapper-utils/compare/v0.4.1...v0.5.0) (2022-08-16)


### Features

* Adding parsing of reference ([#20](https://www.github.com/snakemake/snakemake-wrapper-utils/issues/20)) ([5c259f5](https://www.github.com/snakemake/snakemake-wrapper-utils/commit/5c259f5c9f5d4a15036bc6d2717e23e65cbe5917))
* Regions, samples and targets in get_bcftools_opts ([#19](https://www.github.com/snakemake/snakemake-wrapper-utils/issues/19)) ([e916c8b](https://www.github.com/snakemake/snakemake-wrapper-utils/commit/e916c8b56600655798ff7c1e12133aa44000035c))

### [0.4.1](https://www.github.com/snakemake/snakemake-wrapper-utils/compare/v0.4.0...v0.4.1) (2022-06-30)


### Bug Fixes

* fixed missing return values ([#15](https://www.github.com/snakemake/snakemake-wrapper-utils/issues/15)) ([c6f0fa7](https://www.github.com/snakemake/snakemake-wrapper-utils/commit/c6f0fa71affe7c89a780d44b15050a1ad71cfe1a))
* rename argument in function 'infer_out_format' ([#17](https://www.github.com/snakemake/snakemake-wrapper-utils/issues/17)) ([93d25fa](https://www.github.com/snakemake/snakemake-wrapper-utils/commit/93d25fae4527cbad277aaf9cec4e966b0cc8c81d))

## [0.4.0](https://www.github.com/snakemake/snakemake-wrapper-utils/compare/v0.3.1...v0.4.0) (2022-04-25)


### Features

* function for output format inference and disable temp folder ([#13](https://www.github.com/snakemake/snakemake-wrapper-utils/issues/13)) ([9b88bdd](https://www.github.com/snakemake/snakemake-wrapper-utils/commit/9b88bdd4dccf5a4cecf55facfaca04997b2a3df4))

### [0.3.1](https://www.github.com/snakemake/snakemake-wrapper-utils/compare/v0.3.0...v0.3.1) (2022-01-17)


### Bug Fixes

* migrate publishing to release please ([60a2665](https://www.github.com/snakemake/snakemake-wrapper-utils/commit/60a266593698c5503afbca7e6d5eb21a4c9c0153))
