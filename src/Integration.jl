module Integration

export chebyshev_quadrature, laguerre_quadrature, legendre_quadrature, rectangle, simpsons, trapezoidal

include("Chebyshev.jl")
include("Laguerre.jl")
include("Legendre.jl")
include("Rectangle.jl")
include("Simpson.jl")
include("Trapezoidal.jl")
end # module
