"""
    hermite_quadrature(f::Function, N::Number, k::Integer=1)

when ``k=2`` approximates the integral:

``\\displaystyle \\int_{-\\infty}^{\\infty} e^{-x^2} f(x) dx``

and otherwise approximates the integral:

``\\displaystyle \\int_{-\\infty}^{\\infty} f(x) dx``

in both instances it uses [Gauss-Hermite quadrature](https://en.wikipedia.org/wiki/Chebyshev-Hermite_quadrature) to make this approximation.

``N`` is the number of nodes (or grid points) used. Whilst the type mentioned in the function definition is just 'Number' (as opposed to 'Integer'), this is just so that scientific notation can be used to define it (as scientific notation gives the type 'Float64'), an error message will be printed if it is not a positive integer. 
"""
function hermite_quadrature(f::Function, N::Number, k::Integer=1)
    N = convert(Int64, N);
    nodes, weights = gausshermite(N);
    if k==2
        int = sum(weights.*f.(nodes));
    else
        int = sum(weights.*exp.(nodes.^2).*f.(nodes));
    end
    return int
end
