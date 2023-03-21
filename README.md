# pgmfactors

This project implements discrete factor operations useful for probabilistic graphical models. These are useful in the construction of Hidden Markov Models and Bayesian Networks.

The current version of pgmfactors uses the XTensor library's xt::xarray, an n-dimensional array with shape determined at runtime. This is suited for rapid prototyping, model exploration, and structure learning. Later versions of pgmfactors may include xt::xtensor and xt::xtensor_fixed in order to optimize graphical models that are determined at compile time.

# Building and installing

See the [BUILDING](BUILDING.md) document.

# Licensing

Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or copy at https://www.boost.org/LICENSE_1_0.txt)

