# RVSC-D-26-01958 R reproducibility package — final verified revision

Manuscript: **Pressure-dependent disagreement between oscillometric and invasive arterial pressure measurements in anaesthetized pigs: a retrospective method-comparison study**

This package contains the finalized R reproduction scripts used to verify the revised analyses. The animal-level source data are not distributed because they are corporate records of Boston Scientific Japan K.K.

## Verified environment

- R 4.5.3 (Windows 11 x64)
- readr 2.2.0
- dplyr 1.2.1
- geepack 1.3.13
- lme4 2.1.0
- pROC 1.19.1

## Run order

1. Place authorized data files in `data/`.
2. Run `source("RUN_CORE.R")` for the non-bootstrap analyses.
3. Run `source("RUN_BOOTSTRAP.R")` for the final pig-cluster bootstraps. The LoA bootstrap uses two local workers on Windows.

## Final bootstrap settings

- Deming λ=1: 3,000 pig-cluster resamples, seed 20261104. Verified slope 0.8662667; 95% percentile CI 0.7522364–0.9978518.
- Threshold metrics: 3,000 pig-cluster resamples, seed 20261104.
- Regression-based LoA: 1,000 pig-cluster resamples; bootstrap replicate `b` uses seed `20261104+b`; two workers. At paired mean MAP 120 mmHg, lower LoA point estimate −52.03 mmHg and bootstrap 95% CI approximately −62.2 to −43.8 mmHg.

## GEE ordering

`geepack::geeglm()` requires observations within a cluster to be contiguous. The scripts therefore explicitly sort by pig before pig-clustered GEE models and by session before session-clustered GEE models. The session AR(1) sensitivity uses the validated fast Gaussian GEE implementation matching the sequential-adjacency estimator used when WP3 was frozen.

See `expected/R_FINAL_VERIFICATION_SUMMARY.csv` for the verified target values and `expected/sessionInfo_final.txt` for the frozen software environment.


## Citation and licence

This repository is released under the MIT License. Citation metadata are provided in `CITATION.cff`. For the manuscript, cite the version-specific Zenodo DOI assigned to the archived GitHub release.

## Public repository contents

Only analysis scripts, documentation, verified aggregate results, and the frozen software-environment record are included. Animal-level source data are not included. The `.gitignore` file excludes files placed under `data/` other than `data/README.md` to reduce the risk of accidental source-data publication.
