"""
    chebyshev_quadrature(f::Function, k::Integer, N::Number, a::Number, b::Number)

Uses Chebyshev-Gauss quadrature to approximate the integral:

``\\displaystyle\\int_a^b f(x) dx``.

`N` is the number of nodes (or grid points) used.

`k` is the kind of the Chebyshev polynomial used in the quadrature. If ``k=2``, Chebyshev polynomials of the second kind ``U_n(x)`` are used; if ``k=1``, Chebyshev polynomials of the first kind ``T_n(x)`` are used.
"""
function chebyshev_quadrature(f::Function, k::Integer, N::Number, a::Number, b::Number)
    N = convert(Int64, N);
    n = 1:1:N;
    if k == 2
        # Grid
        x = -cos.(pi*n./(N+1));
        # Linear transformation of the original integration interval [-1,1]
        # to [a,b]
        u = ((b-a)/2)*x.+(a+b)/2;
        # weights are pi/(N+1)*sin^2(pi*n/(N+1)), which equals
        # pi/(N+1)*(1-x^2)
        int = ((b-a)*pi)/(2*(N+1))*sum((sqrt.((-x.^2).+1)).*f.(u));
    elseif k == 1
        # Grid
        x = -cos.(((2*n.-1)*pi)/(2*N));
        # Linear transformation of the original integration interval [-1,1]
        # to [a,b]
        u = ((b-a)/2)*x.+(a+b)/2;
        int = (((b-a)*pi)/(2*N))*sum((sqrt.((-x.^2).+1)).*f.(u));
    end
    return int
end
