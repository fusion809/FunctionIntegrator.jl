"""
    radau_quadrature(f::Function, N::Number, a::Number, b::Number)

uses [Radau quadrature](https://mathworld.wolfram.com/RadauQuadrature.html) to approximate:

``\\displaystyle \\int_a^b f(x) dx.``

``N`` is the number of nodes used in the approximation.
"""
function radau_quadrature(f::Function, N::Number, a::Number, b::Number)
    N = convert(Int64, N);
    zeros, weights = gaussradau(N);
    u = (b-a)/2*zeros.+(a+b)/2;
    int = (b-a)/2*sum(weights.*f.(u));

    return int
end
