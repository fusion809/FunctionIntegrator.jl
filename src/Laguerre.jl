"""
    laguerre_quadrature(f::Function, N::Number, k::Integer=1)

approximates the integral:

``\\displaystyle \\int_0^{\\infty} f(x) dx``

or if k is set to something other than 1:

``\\displaystyle \\int_0^{\\infty} e^{-x} f(x) dx``

using [Gauss-Laguerre quadrature](https://en.wikipedia.org/wiki/Gauss-Laguerre_quadrature).

`N` is the number of nodes (or grid points) used. Whilst the type mentioned in the function definition is just `Number` (as opposed to 'Integer'), this is just so that scientific notation can be used to define it (as scientific notation gives the type 'Float64'), an error message will be printed if it is not a positive integer. 
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
