# Changelog

Version numbers have the form MAJOR.MINOR.PATCH. MAJOR changes when existing calls may stop working.

## 1.0.0 (not released yet)

First numbered version.

- Version function `stm_version`; `stm_dependencies` checks the required versions
- Fixed demo scripts and `es_default_options`, tests and README
- `es_default_options` removed: the identical function in the Metabolic Network Toolbox is used
- `stm_sample_model` and `stm_sample_elasticities` check the dependencies automatically, once per Matlab session
- All functions and folders starting with `es_` (and the folder `es-utils`) renamed to `stm_` (`stm-utils`). This does not apply to `es_default_options`, which is part of the Metabolic Network Toolbox

The state before the cleanup of October 2026 is available under the git tag `v0.9-legacy`.
