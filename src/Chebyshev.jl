"""
    chebyshev_quadrature(f::Function, N::Number, k::Integer, a::Number, b::Number)

uses [Chebyshev-Gauss quadrature](https://en.wikipedia.org/wiki/Chebyshev-Gauss_quadrature) to approximate the integral:

``\\displaystyle\\int_a^b f(x) dx``.

`N` is the number of nodes (or grid points) used.

`k` is the kind of the Chebyshev polynomial used in the quadrature. If ``k=2``, Chebyshev polynomials of the second kind ``U_n(x)`` are used; if ``k=1``, Chebyshev polynomials of the first kind ``T_n(x)`` are used.
"""
function chebyshev_quadrature(f::Function, N::Number, k::Integer, a::Number, b::Number)
    N = convert(Int64, N);
    zeros, weights = gausschebyshev(N, k);
    if k == 2
        # Grid
        # Linear transformation of the original integration interval [-1,1]
        # to [a,b]
        u = ((b-a)/2)*zeros.+(a+b)/2;
        # weights are pi/(N+1)*sin^2(pi*n/(N+1)), which equals
        # pi/(N+1)*(1-x^2)
        int = (b-a)/2*sum(weights.*f.(u).*(sqrt.((-zeros.^2).+1)).^(-1));
    elseif k == 1
        # Linear transformation of the original integration interval [-1,1]
        # to [a,b]
        u = ((b-a)/2)*zeros.+(a+b)/2;
        int = (b-a)/2*sum(weights.*(sqrt.((-zeros.^2).+1)).*f.(u));
    else
        println("k must be either 1 or 2")
    end
    return int
end
