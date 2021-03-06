function stepwise_simpsons(f::Function, h::Number, x::Number)
    return h/6 * (f(x) + 4*f(x+h/2) + f(x+h));
end

function stepwise_simpsons38(f::Function, h::Number, x::Number, i::Integer, N::Number)
    if i==1 || i == N + 1
        return (3*h/8)*f(x)
    elseif ((i-1) % 3) == 0
        return (3*h/4)*f(x)
    else
        return (9*h/8)*f(x)
    end
end
"""
    simpsons_rule(f::Function, N::Number, a::Number, b::Number)

uses [Simpson's rule](https://en.wikipedia.org/wiki/Simpson%27s_rule) to approximate:

``\\displaystyle\\int_a^b f(x) dx``.

``N+1`` steps are used, if endpoints are included. Whilst the type mentioned in the function definition for ``N`` is just 'Number' (as opposed to 'Integer'), this is just so that scientific notation can be used to define it (as scientific notation gives the type 'Float64').
"""
function simpsons_rule(f::Function, N::Number, a::Number, b::Number)
    N = convert(Int64, N);
    h = (b-a)/N;
    y = 0;
    x = a;
    for i=1:N
        y += stepwise_simpsons(f, h, x);
        x += h;
    end
    return y
end

"""
    simpsons38_rule(f::Function, N::Number, a::Number, b::Number)

uses [Simpson's 3/8 rule](https://en.wikipedia.org/wiki/Simpson%27s_rule#Simpson's_3/8_rule) to approximate:

``\\displaystyle\\int_a^b f(x) dx``.

``N+1`` steps are used, if endpoints are included. Whilst the type mentioned in the function definition for ``N`` is just 'Number' (as opposed to 'Integer'), this is just so that scientific notation can be used to define it (as scientific notation gives the type 'Float64'), an error message will be printed if it is not a positive integer. In order to apply Simpson's 3/8 rule, ``N`` must be divisible by 3.
"""
function simpsons38_rule(f::Function, N::Number, a::Number, b::Number)
    N = convert(Int64, N);
    (N % 3) == 0 || error("The input N in simpsons38_rule(f::Function, N::Number, a::Number, b::Number) must be divisible by 3 in order for it to be possible to apply Simpson's 3/8 rule.")
    h = (b-a)/N;
    y = 0;
    x = a;
    for i=1:N+1
        y = y + stepwise_simpsons38(f, h, x, i, N);
        if i < N+1
            x += h;
        end
    end
    return y
end

"""
    adaptive_simpsons_rule(f::Function, a::Number, b::Number, ε::Float64=1e-10)

uses [adaptive Simpson's rule](https://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method) with Lyness's criterion to approximate the integral:

``\\displaystyle \\int_a^b f(x) dx``

with a relative tolerance of ``\\epsilon`` (which by default is 1e-10).
"""
function adaptive_simpsons_rule(f::Function, a::Number, b::Number, ε::Float64=1e-10)
    N = 10;
    criterion = 100*ε;
    while criterion >= 15*ε
        criterion = abs((simpsons_rule(f, N, a, (a+b)/2)+simpsons_rule(f, N, (a+b)/2, b)-simpsons_rule(f, N, a, b))/simpsons_rule(f, N, a, b));
        N = 2*N;
    end
    return simpsons_rule(f, N, a, b)
end
