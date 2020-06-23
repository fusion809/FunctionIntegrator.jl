module Integration

export chebyshev_quadrature, legendre_quadrature, simpsons, trapezoidal

include("Chebyshev.jl")
include("Legendre.jl")
include("Simpson.jl")
include("Trapezoidal.jl")
end # module
