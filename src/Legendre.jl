"""
    legendre(d, n, x)

Computes the [Legendre polynomial](https://mathworld.wolfram.com/LegendrePolynomial.html) of nth order or its derivative, depending on the value of `d`. If ``d=0``, ``P_n(x)`` is returned; if ``d=1``, ``P_n'(x)`` is returned.
"""

function legendre(d, n, x)
    # d is the order of the derivative, if d==0 we're talking about
    # the Legendre poly itself
    # only the 1st order derivative can be found
    if d == 1 && n==0
        return 0
    elseif d==1 && n==1
        return 1
    elseif d==2 && n==0
        return 0
    elseif d==2 && n==1
        return 0
    elseif d==2 && n==2
        return 3
    elseif d==0 && n==0
        return 1
    elseif d==0 && n==1
        return x
    else
        # Check if x is a scalar and define Pn or Pn' accordingly
        if typeof(x) == Float64 || typeof(x) == Int64
            Pn = ones(n+1);
            Pn[2] = x;
            for i=2:n
                k = i-1;
                Pn[i+1] = ((2*k+1)*x*Pn[i]-k*Pn[i-1])/(k+1);
            end
            if d==0
                return Pn[n+1]
            elseif d==1
                return n*(x*Pn[n+1]-Pn[n])./((x.^2).-1)
            elseif d==2
                return ((x^2*((n+1)*(n+2)-2*(2*n+1))-n*(n+1))*Pn[n+1]+2*n*x.*Pn[n]).*((x.^2).-1).^(-2)
            end
        # x still needs to be a vector of some description.
        else
            N = length(x);
            Pn = ones(N,n+1);
            Pn[:,2] = x;
            for i=2:n
                k = i-1;
                Pn[:,i+1] = ((2*k+1)*x.*Pn[:,i].-k*Pn[:,i-1])/(k+1);
            end
            if d==0
                return Pn[:,n+1]
            elseif d==1
                return n*(x.*Pn[:,n+1].-Pn[:,n])./((x.^2).-1)
            elseif d==2
                return ((x^2*((n+1)*(n+2)-2*(2*n+1))-n*(n+1))*Pn[:,n+1]+2*n*x.*Pn[:,n]).*((x.^2).-1).^(-2)
            end
        end
    end
end

"""
    legendre_zeros(n)

Returns the zeros of the nth [Legendre polynomial](https://mathworld.wolfram.com/LegendrePolynomial.html).
"""
function legendre_zeros(n)
    N = 2*n;
    # Chebyshev roots grid, seems to be the best analytical approximator
    # of the Legendre roots
    x=-cos.(pi*(2*(1:N).-1)/(2*N));
    xx=ones(length(x));
    roots = ones(n);
    k = 1;

    for i in 2:length(x)
        # xx is where our zeros will initially go
        xx[i] = x[i];
        diff = -legendre(0,n,xx[i])/legendre(1,n,xx[i]);

        # Use Newton's method until diff </= 1e-12
        while abs(diff) > 1e-12
            diff = -legendre(0,n,xx[i])/legendre(1,n,xx[i]);
            xx[i] += diff;
        end

        # The next bit of code is necessary to ensure there is no doubling
        # up in our roots
        Diff = 1e-5;
        if i > 1
            for j in 1:i-1
                # Check xx differs from other xx values, hence the root
                # is unique
                if abs(xx[i]-xx[j]) < Diff
                    Diff = abs(xx[i]-xx[j]);
                end
            end
            # Some elements of xx will be NaN, so let's filter those out
            # and those elements that are identical
            if Diff > 1e-10 && !isnan(xx[i])
                roots[k] = xx[i];
                k += 1;
            end
        end
    end
    # roots is where our unique zeros go
    return roots = sort(roots)
end

"""
    weights, zeros = legendre_weights_and_zeros(n)

Returns the weights and nodes for Legendre-Gauss quadrature with ``n`` nodes.
"""

# weights, zeros = legendre_weights_and_zeros(n)
# where zeros are the zeros of the nth order Legendre polynomial
# and weights are the weights for nth order Legendre quadrature
function legendre_weights_and_zeros(n)
    zeros = legendre_zeros(n);
    # Cannot use 2/(((-zeros.^2).+1).*L.^2) as otherwise it seems to use
    # least squares division
    return 2*((((-zeros.^2).+1).*legendre(1,n,zeros).^2).^(-1)), zeros
end

"""
    legendre_quadrature(f, N, a, b)

Uses [Legendre-Gauss quadrature](https://mathworld.wolfram.com/Legendre-GaussQuadrature.html) to approximate:

``\\displaystyle \\int_a^b f(x) dx``.

The argument `N` refers to the number of nodes used in the quadrature.
"""
function legendre_quadrature(f, N, a, b)
    N = convert(Int64, N);
    weights, zeros = legendre_weights_and_zeros(N);
    u = (b-a)/2*zeros.+(a+b)/2;
    return (b-a)/2*sum(weights.*f(u))
end
