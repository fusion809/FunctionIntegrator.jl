# FunctionIntegrator.jl
![Travis](https://travis-ci.com/fusion809/FunctionIntegrator.jl.svg?branch=master)

This package should be treated as a second-rate alternative to the excellent [QuadGK](https://github.com/JuliaMath/QuadGK.jl) package. QuadGK provides more accurate integration for many problems, and also provides an error estimate which functions in this package do not.

This package provides the following functions:

* `chebyshev_quadrature(f, N, k, a, b)`
* `hermite_quadrature(f, N, k)`
* `jacobi_quadrature(f, N, α, β, a, b)`
* `laguerre_quadrature(f, N, k)`
* `legendre_quadrature(f, N, a, b)`
* `lobatto_quadrature(f, N, a, b)`
* `radau_quadrature(f, N, a, b)`
* `rectangle(f, N, a, b)`
* `simpsons(f, N, a, b)`
* `trapezoidal(f, N, a, b)`

use Julia's help function to find out usage information, should you need it.

This package is not currently in Julia's official repositories, so to install it one has to run:

```julia
pkg> add git://github.com/fusion809/Integration.jl.git
```

## Choice of function
As a general rule of thumb, `simpsons` should be the function you use when you're unsure which function to use as its approximations with large N (e.g. 1e8) are nearly always accurate to at least 4 digits. Despite this, for many problems some of the `_quadrature` functions may provide more accurate results with a far smaller N value. The main time when `simpsons` should be avoided is when there are singularities at the endpoints of the domain of integration, in which case using `chebyshev_quadrature` with k=1 or using `legendre_quadrature` is likely best.

Each of the functions whose name ends with `_quadrature` use [Gaussian quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature), the specifics of which differ between functions. Out of them, `chebyshev_quadrature` and `legendre_quadrature` are perhaps the best to go with when you're uncertain which out of the `_quadrature` functions to go with as they can be applied to any grid and any problem and *usually* arrive at a result that is accurate to at least 4 decimal places provided N&geq;1e6. These two functions have the added benefit of being faster than all other `_quadrature` functions, after controlling for N.

## Acknowledgements
I'd like to thank the Julia discourse community for their generous help through my Julia journey, and I would also like to thank the developers of  Julia, as without their work this package would not even be possible (I know, obviously, as this is a Julia package). Likewise, I'd also like to thank the developers of [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl), as without their efficient algorithms for finding the nodes and weights for various Gaussian quadrature techniques, several of the functions in this repository would be far less efficient (especially `legendre_quadrature`), or may not even exist. I'd also like to thank the developers of [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl), as some of their functions are useful for testing the functions in this package.
