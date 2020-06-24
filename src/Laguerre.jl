using FastGaussQuadrature

"""
    laguerre_quadrature(f::Function, N::Number, k::Integer=1)

approximate the integral:

``\\displaystyle \\int_0^{\\infty} f(x) dx``

or if k is set to something other than 1:

``\\displaystyle \\int_0^{\\infty} e^{-x} f(x) dx``

using Gauss-Laguerre quadrature. ``N`` is the number of quadrature nodes used.
"""
function laguerre_quadrature(f::Function, N::Number, k::Integer=1)
    N = convert(Int64, N);
    zeros, weights = gausslaguerre(N);
    if k==1
        int = sum(exp.(zeros).*f.(zeros).*weights);
    else
        int = sum(f.(zeros).*weights);
    end
    return int;
end
