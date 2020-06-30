"""
    jacobi_quadrature(f::Function, N::Number, α::Number, β::Number, a::Number, b::Number)

uses [Gauss-Jacobi quadrature](https://en.wikipedia.org/wiki/Gauss-Jacobi_quadrature) to approximate:

``\\displaystyle \\int_a^b f(x) dx.``

The remaining inputs of this function are:
* ``N`` is the number of nodes (or grid points) used. Whilst the type mentioned in the function definition is just 'Number' (as opposed to 'Integer'), this is just so that scientific notation can be used to define it (as scientific notation gives the type 'Float64'), an error message will be printed if it is not a positive integer. 
* ``α`` and ``β``, which are parameters of the weighting function.
"""
function jacobi_quadrature(f::Function, N::Number, α::Number, β::Number, a::Number, b::Number)
    N              = convert(Int64, N);
    nodes, weights = gaussjacobi(N, α, β);
    u              = (b-a)/2*nodes.+(a+b)/2;
    int            = (b-a)/2*sum(weights.*f.(u).*(-nodes.+1.0).^(-α).*(nodes.+1.0).^(-β));

    return int
end
