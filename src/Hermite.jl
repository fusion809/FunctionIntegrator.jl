"""
    hermite_quadrature(f::Function, N::Number, k::Integer=1)

when ``k=2`` approximates the integral:

``\\displaystyle \\int_{-\\infty}^{\\infty} e^{-x^2} f(x) dx``

and otherwise approximates the integral:

``\\displaystyle \\int_{-\\infty}^{\\infty} f(x) dx``

in both instances it uses [Gauss-Hermite quadrature](https://en.wikipedia.org/wiki/Chebyshev-Hermite_quadrature) to make this approximation. ``N`` quadrature nodes are used in this approximation.
"""
function hermite_quadrature(f::Function, N::Number, k::Integer=1)
    N = convert(Int64, N);
    zeros, weights = gausshermite(N);
    if k==2
        int = sum(weights.*f.(zeros));
    else
        int = sum(weights.*exp.(zeros.^2).*f.(zeros));
    end
    return int
end
