function stepwise_trapezoidal(f, h, x, i, N)
    return h/2 * (f(x) + f(x+h));
end

"""
    trapezoidal_rule(f::Function, N::Number, a::Number, b::Number)

uses the [trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule) to approximate:

``\\displaystyle\\int_a^b f(x) dx``.

``N+1`` steps are used, if endpoints are included. Whilst the type mentioned in the function definition for ``N`` is just 'Number' (as opposed to 'Integer'), this is just so that scientific notation can be used to define it (as scientific notation gives the type 'Float64'), an error message will be printed if it is not a positive integer. 
"""
function trapezoidal_rule(f::Function, N::Number, a::Number, b::Number)
    N = convert(Int64, N);
    h = (b-a)/N;
    y = 0;
    x = a;
    for i=1:N
        y += stepwise_trapezoidal(f, h, x, i, N);
        x += h;
    end
    return y
end
