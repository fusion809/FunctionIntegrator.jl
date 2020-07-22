"""
    rombergs_method(f::Function, N::Number, a::Number, b::Number)

computes the integral:

``\\displaystyle \\int_a^b f(x) dx``

using [Romberg's method](https://en.wikipedia.org/wiki/Romberg's_method#Method) with the ``n`` mentioned in the linked article being equal to ``N`` in the function arguments. 
"""
function rombergs_method(f::Function, N::Number, a::Number, b::Number)
    N = convert(Int64, N);
    n = 1:N;
    h = (2.0 .^ (-n)).*(b-a);
    n = nothing;
    R = zeros(N+1,N+1);
    # This is really R[0,0]
    R[1,1] = h[1]*(f.(a)+f.(b));
    for n=2:N+1
        upper_bounds = convert(Int64, 2^(n-2));
        k = 1:upper_bounds;
        R[n,1] = 1/2*R[n-1,1] + h[n-1]*sum(f.((2*k.-1)*h[n-1].+a));
        k = nothing;
        for m=2:n
            R[n,m] = (4.0^(m-1) - 1).^(-1)*((4^(m-1)*R[n,m-1]-R[n-1,m-1]));
        end
    end
    return R[N+1,N+1]
end