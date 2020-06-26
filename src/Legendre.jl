"""
    legendre_quadrature(f::Function, N::Number, a::Number, b::Number)

Uses [Legendre-Gauss quadrature](https://en.wikipedia.org/wiki/Gauss-Legendre_quadrature) to approximate:

``\\displaystyle \\int_a^b f(x) dx``.

``N`` is the number of nodes (or grid points) used. Whilst the type mentioned in the function definition is just 'Number' (as opposed to 'Integer'), this is just so that scientific notation can be used to define it (as scientific notation gives the type 'Float64'), an error message will be printed if it is not a positive integer. 
"""
function legendre_quadrature(f::Function, N::Number, a::Number, b::Number)
    N = convert(Int64, N);
    zeros, weights = gausslegendre(N);
    u = (b-a)/2*zeros.+(a+b)/2;
    return (b-a)/2*sum(weights.*f.(u))
end
