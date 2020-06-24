module Integration

export chebyshev_quadrature, legendre_quadrature, rectangle, simpsons, trapezoidal

include("Chebyshev.jl")
include("Legendre.jl")
include("Rectangle.jl")
include("Simpson.jl")
include("Trapezoidal.jl")
end # module
