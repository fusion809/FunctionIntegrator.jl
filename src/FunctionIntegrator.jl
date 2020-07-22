module FunctionIntegrator

using FastGaussQuadrature

export chebyshev_quadrature, hermite_quadrature, jacobi_quadrature, laguerre_quadrature, legendre_quadrature, lobatto_quadrature, radau_quadrature, rectangle_rule_left, rectangle_rule_midpoint, rectangle_rule_right, rombergs_method, simpsons_rule, simpsons38_rule, adaptive_simpsons_rule, trapezoidal_rule

include("Chebyshev.jl")
include("Hermite.jl")
include("Jacobi.jl")
include("Laguerre.jl")
include("Legendre.jl")
include("Lobatto.jl")
include("Radau.jl")
include("Rectangle.jl")
include("Romberg.jl")
include("Simpson.jl")
include("Trapezoidal.jl")
end # module
