# Integration.jl

This package should be treated as a second-rate alternative to the excellent [QuadGK](https://github.com/JuliaMath/QuadGK.jl) package. QuadGK provides more accurate integration for many problems, and also provides an error estimate which functions in this package do not.

This package provides the following functions:

* `chebyshev_quadrature(f, k, N, a, b)`
* `legendre_quadrature(f, N, a, b)`
* `simpsons(f, N, a, b)`
* `trapezoidal(f, N, a, b)`

use Julia's help function to find out usage information, should you need it. Out of these, I personally recommend you use `simpsons` when you're unsure which to use, as it nearly always gives an approximately accurate result, whilst the two `_quadrature` functions sometimes do not (in fact, they sometimes converge on an answer that is incredibly inaccurate, see [test/runtests.jl#L59](test/runtests.jl#L59) for an example function for which this is the case).
