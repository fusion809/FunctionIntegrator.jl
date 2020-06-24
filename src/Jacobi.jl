"""
    jacobi_quadrature(f::Function, N::Number, α::Number, β::Number, a::Number, b::Number)

uses Gauss-Jacobi quadrature to approximate:

``\\int_a^b f(x) dx.``

The remaining inputs of this function are:
* ``N``, which refers to the number of nodes used in the quadrature.
* ``α`` and ``β``, which are parameters of the weighting function.
"""
function jacobi_quadrature(f::Function, N::Number, α::Number, β::Number, a::Number, b::Number)
    N              = convert(Int64, N);
    zeros, weights = gaussjacobi(N, α, β);
    u              = (b-a)/2*zeros.+(a+b)/2;
    int            = (b-a)/2*sum(weights.*f.(u).*(-zeros.+1.0).^(-α).*(zeros.+1.0).^(-β));

    return int
end
