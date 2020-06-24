# Integration.jl
![Travis](https://travis-ci.com/fusion809/Integration.jl.svg?branch=master)

This package should be treated as a second-rate alternative to the excellent [QuadGK](https://github.com/JuliaMath/QuadGK.jl) package. QuadGK provides more accurate integration for many problems, and also provides an error estimate which functions in this package do not.

This package provides the following functions:

* `chebyshev_quadrature(f, k, N, a, b)`
* `legendre_quadrature(f, N, a, b)`
* `rectangle(f, N, a, b)`
* `simpsons(f, N, a, b)`
* `trapezoidal(f, N, a, b)`

use Julia's help function to find out usage information, should you need it.
