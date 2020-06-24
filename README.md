# Integration.jl
![Travis](https://travis-ci.com/fusion809/Integration.jl.svg?branch=master)

This package should be treated as a second-rate alternative to the excellent [QuadGK](https://github.com/JuliaMath/QuadGK.jl) package. QuadGK provides more accurate integration for many problems, and also provides an error estimate which functions in this package do not.

This package provides the following functions:

* `chebyshev_quadrature(f, k, N, a, b)`
* `hermite_quadrature(f, N, k)`
* `laguerre_quadrature(f, N, k)`
* `legendre_quadrature(f, N, a, b)`
* `rectangle(f, N, a, b)`
* `simpsons(f, N, a, b)`
* `trapezoidal(f, N, a, b)`

use Julia's help function to find out usage information, should you need it.

This package is not currently in Julia's official repositories, so to install it one has to run:

```julia
pkg> add git://github.com/fusion809/Integration.jl.git
```

## Acknowledgements
I'd like to thank the Julia discourse community for their generous help through my Julia journey, and I would also like to thank the developers of  Julia, as without their work this package would not be possible. Likewise, I'd also like to thank the developers of [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl), as without their efficient algorithms for finding the nodes and weights for various Gaussian quadrature techniques, several of the functions in this repository would be far less efficient (especially `legendre_quadrature`), or may not even exist. I'd also like to thank the developers of [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl), as some of their functions are useful for testing the functions in this package.
