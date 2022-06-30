  <!-- badges: start -->
[![R-CMD-check](https://github.com/adamsardar/metaDEGth/workflows/R-CMD-check/badge.svg)](https://github.com/adamsardar/metaDEGth/actions)
[![Codecov test coverage](https://codecov.io/gh/adamsardar/metaDEGth/branch/master/graph/badge.svg)](https://app.codecov.io/gh/adamsardar/metaDEGth?branch=master)
[![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-blue.svg)](https://www.tidyverse.org/lifecycle/#experimental)
  <!-- badges: end -->

Following the approach of Morris & Pounds, a distribution of P-values can be decomposed as a mixture of two β-distributions, one with shape paramter (α) of 1 - modelling noise under the null-hypothesis - and another with an emperically fit shape parament (α) of a,  modelling the signal component. One can view this as modeling P-values as a random mixture process, where P ~ (1-λ)β(a,1) + λβ(1,1).

[Fisher's method](https://en.wikipedia.org/wiki/Fisher%27s_method) of combining independent P-values can be used to conduct a meta-analysis across several tests. However, it assumes a null hypothesis of uniformly distributed P-values, which in the presence of signal is manifestly incorrect. This work looks to produce an analagous test against an appropriately specified null-hypothesis that respects the pre-existing signal.

The aim of this work is a short paper for BioArXiv/PeerJ and a small R package.

Yes, the package name is a bad play on the [80's metal band.](https://en.wikipedia.org/wiki/Megadeth)