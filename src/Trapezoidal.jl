function stepwise_trapezoidal(f, h, x, i, N)
    if i == 1 || i == N + 1
        return h / 2 * f(x)
    else
        return h * f(x)
    end
end

"""
    trapezoidal(f, N, a, b)

Uses the [trapezoidal rule](https://mathworld.wolfram.com/TrapezoidalRule.html) to approximate:

``\\displaystyle\\int_a^b f(x) dx``.

``N+1`` steps are used, if endpoints are included. 
"""
function trapezoidal(f, N, a, b)
    h = (b-a)/N;
    y = 0;
    x = a;
    for i=1:N+1
        y = y + stepwise_trapezoidal(f, h, x, i, N);
        if i < N+1
            x = x + h;
        end
    end
    return y
end
