# Changelog

## [0.8.0](https://www.github.com/snakemake/snakemake-wrapper-utils/compare/v0.7.2...v0.8.0) (2025-09-11)


### Features

* Mem overhead ([#47](https://www.github.com/snakemake/snakemake-wrapper-utils/issues/47)) ([fdce4cd](https://www.github.com/snakemake/snakemake-wrapper-utils/commit/fdce4cd716cabe625f99006da86b890aeb407bd5))
* move temp output to final destinations ([#51](https://www.github.com/snakemake/snakemake-wrapper-utils/issues/51)) ([fd1df7e](https://www.github.com/snakemake/snakemake-wrapper-utils/commit/fd1df7e3a43f70bb15e855e70656e23819029304))


### Bug Fixes

* Add more FASTA extensions ([#49](https://www.github.com/snakemake/snakemake-wrapper-utils/issues/49)) ([2cce124](https://www.github.com/snakemake/snakemake-wrapper-utils/commit/2cce1247feee4a1468f138abfde6f86c04f3a2ff))

### [0.7.2](https://www.github.com/snakemake/snakemake-wrapper-utils/compare/v0.7.1...v0.7.2) (2025-03-21)


### Bug Fixes

* delete missing parenthesis ([#45](https://www.github.com/snakemake/snakemake-wrapper-utils/issues/45)) ([5d3bbb1](https://www.github.com/snakemake/snakemake-wrapper-utils/commit/5d3bbb155db575eb8818a399046089b541738de7))

### [0.7.1](https://www.github.com/snakemake/snakemake-wrapper-utils/compare/v0.7.0...v0.7.1) (2025-03-21)


### Bug Fixes

* typo on samtools regions file ([#43](https://www.github.com/snakemake/snakemake-wrapper-utils/issues/43)) ([790f5fa](https://www.github.com/snakemake/snakemake-wrapper-utils/commit/790f5faace2eaa7541d122170c2a4561fb88c810))

## [0.7.0](https://www.github.com/snakemake/snakemake-wrapper-utils/compare/v0.6.2...v0.7.0) (2025-03-10)


### Features

* GATK utils for common IO files ([#39](https://www.github.com/snakemake/snakemake-wrapper-utils/issues/39)) ([4ba82aa](https://www.github.com/snakemake/snakemake-wrapper-utils/commit/4ba82aa2ea3c2f675f94b9d1a09f7ab0b33c3155))
* parse file format from extension (ignore compression) ([#41](https://www.github.com/snakemake/snakemake-wrapper-utils/issues/41)) ([b0820f2](https://www.github.com/snakemake/snakemake-wrapper-utils/commit/b0820f24169d457fa2b1b7a16168b509dc787720))
* parse samtools regions ([#40](https://www.github.com/snakemake/snakemake-wrapper-utils/issues/40)) ([2fa70b8](https://www.github.com/snakemake/snakemake-wrapper-utils/commit/2fa70b89b529e3d0b0a1fee38c4479b634342b68))

### [0.6.2](https://www.github.com/snakemake/snakemake-wrapper-utils/compare/v0.6.1...v0.6.2) (2023-07-06)


### Bug Fixes

* bug where mem_mb is never read ([#34](https://www.github.com/snakemake/snakemake-wrapper-utils/issues/34)) ([be7b4b2](https://www.github.com/snakemake/snakemake-wrapper-utils/commit/be7b4b2f6ad4afb88fb01aae4152ec21967f32bf))

### [0.6.1](https://www.github.com/snakemake/snakemake-wrapper-utils/compare/v0.6.0...v0.6.1) (2023-06-29)


### Bug Fixes

* set default memory for `get_mem()` to 1/5 of 1GB ([#32](https://www.github.com/snakemake/snakemake-wrapper-utils/issues/32)) ([999240e](https://www.github.com/snakemake/snakemake-wrapper-utils/commit/999240ebb27adf89ed4565e54c4171a2885123ee))

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
