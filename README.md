# FunctionIntegrator.jl
![Travis](https://travis-ci.com/fusion809/Integration.jl.svg?branch=master)
![CompatHelper](https://github.com/fusion809/FunctionIntegrator.jl/workflows/CompatHelper/badge.svg?event=push)

This package should be treated as a second-rate alternative to the excellent [QuadGK](https://github.com/JuliaMath/QuadGK.jl) package. QuadGK provides more accurate integration for many problems, and also provides an error estimate which functions in this package do not.

This package provides the following functions:

* `chebyshev_quadrature(f, N, k, a, b)`
* `hermite_quadrature(f, N, k)`
* `jacobi_quadrature(f, N, α, β, a, b)`
* `laguerre_quadrature(f, N, k)`
* `legendre_quadrature(f, N, a, b)`
* `lobatto_quadrature(f, N, a, b)`
* `radau_quadrature(f, N, a, b)`
* `rectangle_rule(f, N, a, b)`
* `simpsons_rule(f, N, a, b)`
* `trapezoidal_rule(f, N, a, b)`

use Julia's help function (e.g. by typing `?chebyshev_quadrature`) to find out usage information, should you need it.

This package is not currently in Julia's registries, so to install it one has to run:

```julia
pkg> add git://github.com/fusion809/FunctionIntegrator.jl.git
```

Although, a [pull request](https://github.com/JuliaRegistries/General/pull/16929) to add this package to the "General" registry is currently open, so hopefully soon this will no longer be necessary.

## Choice of function
As a general rule of thumb, `simpsons_rule` should be the function you use when you're unsure which function to use as its approximations with large N (e.g. 1e8) are nearly always accurate to at least 6 digits. Despite this, for many problems some of the `_quadrature` functions may provide more accurate results with a far smaller N value. The main time when `simpsons_rule` should be avoided is when there are singularities at the endpoints of the domain of integration, in which case using `chebyshev_quadrature` with k=1 or using `legendre_quadrature` is likely best.

Each of the functions whose name ends with `_quadrature` uses [Gaussian quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature), the specifics of which differ between functions. Out of them, `chebyshev_quadrature` and `legendre_quadrature` are perhaps the best to go with when you're uncertain which out of the `_quadrature` functions to go with as they can be applied to any grid and any problem and *usually* arrive at a result that is accurate to at least 4 decimal places provided N&geq;1e6. These two functions have the added benefit of being faster than all other `_quadrature` functions, after controlling for N.

The [test/](test/) folder has test scripts that approximate various different integrals (each file, except [test/runtests.jl](test/runtests.jl), pertains to a different integral); the N values given are the smallest possible to pass each of the tests listed (except when the test involves a less than (<) sign). If you want to know which function to use for which integral, these tests may be useful as a rough guide.

| Function               | Domain of integration | Arguments | Notes |
|------------------------|-----------------------|-----------|-------|
| `chebyshev_quadrature` | [a,b]                 | `f`, the function being integrated.<br/>`N`, the number of grid points.<br/>`k`, the type of Chebyshev quadrature being used. 1, 2, 3, and 4 refer to the Chebyshev T, U, V and W polynomials respectively.<br/>`a`, the start of the domain of integration.<br/>`b`, the end of the domain of integration. | Out of `_quadrature` functions, it or `legendre_quadrature` are favoured for most problems. |
| `hermite_quadrature`   | [-∞,∞]                | `f`, the function being integrated.<br/>`N`, the number of grid points.<br/>`k` determines the problem being solved (whether e^(-x^2) is assumed to be part of the integrand (k=2) or not). | Only use this if your integration domain is [-∞,∞] and your integrand converges rapidly. |
| `jacobi_quadrature`    | [a,b]                 | `f`, the function being integrated.<br/>`N`, the number of grid points.<br/>`α` and `β` are parameters of the weighting function.<br/>`a`, the start of the domain of integration.<br/>`b`, the end of the domain of integration. | Potentially useful when the integrand includes powers of (1-x) and (1+x). |
| `laguerre_quadrature`  | [0,∞]                 | `f`, the function being integrated.<br/>`N`, the number of grid points.<br/>`k` determines the problem being solved (whether e^(-x) is assumed to be part of the integrand (k=2) or not). | Only use this if your integration domain is [0,∞] and your integrand converges rapidly. |
| `legendre_quadrature`  | [a,b]                 | `f`, the function being integrated.<br/>`N`, the number of grid points.<br/>`a`, the start of the domain of integration.<br/>`b`, the end of the domain of integration. | Generally, this is the best `_quadrature` function to go with, when you're otherwise unsure which to go with. |
| `lobatto_quadrature`   | [a,b]                 | `f`, the function being integrated.<br/>`N`, the number of grid points.<br/>`a`, the start of the domain of integration.<br/>`b`, the end of the domain of integration. | If there are singularities at either or both of the endpoints this method will fail to give an accurate result. |
| `radau_quadrature`     | [a,b]                 | `f`, the function being integrated.<br/>`N`, the number of grid points.<br/>`a`, the start of the domain of integration.<br/>`b`, the end of the domain of integration. | If there are singularities at either or both of the endpoints this method will fail to give an accurate result. |
| `rectangle_rule`       | [a,b]                 | `f`, the function being integrated.<br/>`N`, the number of grid points.<br/>`a`, the start of the domain of integration.<br/>`b`, the end of the domain of integration. | Usually this is the least accurate method. |
| `simpsons_rule`        | [a,b]                 | `f`, the function being integrated.<br/>`N`, the number of grid points.<br/>`a`, the start of the domain of integration.<br/>`b`, the end of the domain of integration. | One of the best functions to use when you're unsure which to use. |
| `trapezoidal_rule`     | [a,b]                 | `f`, the function being integrated.<br/>`N`, the number of grid points.<br/>`a`, the start of the domain of integration.<br/>`b`, the end of the domain of integration. | Usually this is the second-least accurate method. |

## Acknowledgements
I'd like to thank the Julia discourse community for their generous help through my Julia journey, and I would also like to thank the developers of  Julia, as without their work this package would not even be possible (I know, obviously, as this is a Julia package). Likewise, I'd also like to thank the developers of [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl), as without their efficient algorithms for finding the nodes and weights for various Gaussian quadrature techniques, several of the functions in this repository would be far less efficient (especially `legendre_quadrature`), or may not even exist. I'd also like to thank the developers of [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl), as some of their functions are useful for testing the functions in this package.
