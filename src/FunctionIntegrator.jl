module FunctionIntegrator

using FastGaussQuadrature

export chebyshev_quadrature, hermite_quadrature, jacobi_quadrature, laguerre_quadrature, legendre_quadrature, lobatto_quadrature, radau_quadrature, rectangle_rule, simpsons_rule, trapezoidal_rule

include("Chebyshev.jl")
include("Hermite.jl")
include("Jacobi.jl")
include("Laguerre.jl")
include("Legendre.jl")
include("Lobatto.jl")
include("Radau.jl")
include("Rectangle.jl")
include("Simpson.jl")
include("Trapezoidal.jl")
end # module
