using FastGaussQuadrature

"""
    legendre_quadrature(f, N, a, b)

Uses [Legendre-Gauss quadrature](https://mathworld.wolfram.com/Legendre-GaussQuadrature.html) to approximate:

``\\displaystyle \\int_a^b f(x) dx``.

The argument `N` refers to the number of nodes used in the quadrature.
"""
function legendre_quadrature(f, N, a, b)
    N = convert(Int64, N);
    zeros, weights = gausslegendre(N);
    u = (b-a)/2*zeros.+(a+b)/2;
    return (b-a)/2*sum(weights.*f(u))
end
