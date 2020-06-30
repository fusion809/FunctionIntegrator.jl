"""
    rectangle_rule_left(f::Function, N::Number, a::Number, b::Number)

numerically approximates:

``\\displaystyle\\int_a^b f(x) dx``

using the [left rectangle rule](https://en.wikipedia.org/wiki/Riemann_sum#Left_Riemann_sum) with ``N`` gridpoint values. Whilst the type mentioned in the function definition for ``N`` is just 'Number' (as opposed to 'Integer'), this is just so that scientific notation can be used to define it (as scientific notation gives the type 'Float64'); an error message will be printed if it is not a positive integer.
"""
function rectangle_rule_left(f::Function, N::Number, a::Number, b::Number)
    N = convert(Int64, N);
    h = (b-a)/N;
    y = 0;
    x = a;
    for i=1:N
        y += h*f(x)
        x = x + h;
    end
    return y
end

"""
    rectangle_rule_midpoint(f::Function, N::Number, a::Number, b::Number)

numerically approximates:

``\\displaystyle\\int_a^b f(x) dx``

using the [midpoint rectangle rule](https://en.wikipedia.org/wiki/Riemann_sum#Midpoint_rule) with ``N`` gridpoint values. Whilst the type mentioned in the function definition for ``N`` is just 'Number' (as opposed to 'Integer'), this is just so that scientific notation can be used to define it (as scientific notation gives the type 'Float64'); an error message will be printed if it is not a positive integer.
"""
function rectangle_rule_midpoint(f::Function, N::Number, a::Number, b::Number)
    N = convert(Int64, N);
    h = (b-a)/N;
    y = 0;
    x = a+h/2;
    for i=1:N
        y += h*f(x)
        x = x + h;
    end
    return y
end

"""
    rectangle_rule_right(f::Function, N::Number, a::Number, b::Number)

numerically approximates:

``\\displaystyle\\int_a^b f(x) dx``

using the [right rectangle rule](https://en.wikipedia.org/wiki/Riemann_sum#Right_Riemann_sum) with ``N`` gridpoint values. Whilst the type mentioned in the function definition for ``N`` is just 'Number' (as opposed to 'Integer'), this is just so that scientific notation can be used to define it (as scientific notation gives the type 'Float64'); an error message will be printed if it is not a positive integer.
"""
function rectangle_rule_right(f::Function, N::Number, a::Number, b::Number)
    N = convert(Int64, N);
    h = (b-a)/N;
    y = 0;
    x = a+h;
    for i=1:N
        y += h*f(x)
        x = x + h;
    end
    return y
end
