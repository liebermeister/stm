# Changelog

Version numbers have the form MAJOR.MINOR.PATCH. MAJOR changes when existing calls may stop working.

## 1.0.0 (2026-10-31)

First numbered version.

- HTML documentation removed from the repository (it can be regenerated with m2html)
- Version function `stm_version`; `stm_dependencies` checks the required versions
- Fixed demo scripts and `stm_default_options`, tests and README
- `stm_default_options` removed: the identical function in the Metabolic Network Toolbox is used
- `stm_sample_model` and `stm_sample_elasticities` check the dependencies automatically, once per Matlab session
- All functions and folders starting with `es_` (and the folder `es-utils`) renamed to `stm_` (`stm-utils`); the structs `es_options` and `es_constraints` are now `stm_options` and `stm_constraints`

The state before the cleanup of October 2026 is available under the git tag `v0.9-legacy`.
