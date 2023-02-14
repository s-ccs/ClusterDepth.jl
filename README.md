# ClusterDepth


[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://s-ccs.github.io/ClusterDepth.jl/dev/)
[![Build Status](https://github.com/s-ccs/ClusterDepth.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/s-ccs/ClusterDepth.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/s-ccs/ClusterDepth.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/s-ccs/ClusterDepth.jl)

Fast implementation of the ClusterDepth multiple comparison algorithm from Frossard and Renaud [Neuroimage 2022](https://doi.org/10.1016/j.neuroimage.2021.118824)

This is especially interesting to EEG signals. Currently only acts on a single channel/timeseries. Multichannel as discussed in the paper is the next step.

## Quickstart
```julia
using ClusterDepth
pval_corrected = clusterdepth(erpMatrix; τ=2.3,nperm=5000)
```
![grafik](https://user-images.githubusercontent.com/10183650/218683929-5fc27ca0-8076-479e-b359-a212bda0b346.png)

## FWER check
We checked FWER for `troendle(...)` and `clusterdepth(...)` [(link to docs)](https://www.s-ccs.de/ClusterDepth.jl/dev/reference/type1/)

For clustedepth we used 5000 repetitions, 5000 permutations, 200 tests.
|simulation|noise|uncorrected|type|
|---|---|---|---|
|clusterdepth|white|1.0|0.0554|
|clusterdepth|red*|XX|XX|
|troendle|white|XX|XX|
|troendle|red*|XX|XX|

Uncorrected should be 1 - it is very improbable, than none of the 200 tests in one repetition is not significant (we expect 5% to be).

Corrected should be 0.05 (CI-95 [0.043,0.0564])


\* red noise introduces strong correlations between individual trials, thus makes the tests correlated while following the H0.

## Citing
Algorithm published in https://doi.org/10.1016/j.neuroimage.2021.118824 - Frossard & Renaud 2022, Neuroimage

Some functions are inspired by [R::permuco](https://cran.r-project.org/web/packages/permuco/index.html), written by Jaromil Frossard. Note: Permuco is GPL licensed, but Jaromil Frossard released the relevant clusteterdepth functions to me under MIT. Thereby, this repository can be licensed under MIT.

XXX ZENODO LINK
