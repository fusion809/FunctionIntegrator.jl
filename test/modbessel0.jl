# test an unusual trig integral
x = 1;
function dmodified_bessel_1(theta)
    # ^(-1) is used instead of / as / produces a matrix instead of
    # a vector
    2*x/pi*cosh(x*sin(theta))*(cos(theta))^2
end

@testset "modbessel0" begin
    @test chebyshev_quadrature(dmodified_bessel_1, 4942, 1, 0, pi/2) ≈ besseli(x,1)
    @test chebyshev_quadrature(dmodified_bessel_1, 6987, 2, 0, pi/2) ≈ besseli(x,1)
    @test chebyshev_quadrature(dmodified_bessel_1, 6988, 3, 0, pi/2) ≈ besseli(x,1)
    @test chebyshev_quadrature(dmodified_bessel_1, 4941, 4, 0, pi/2) ≈ besseli(x,1)
    @test jacobi_quadrature(dmodified_bessel_1, 10896, 1, 1, 0, pi/2) ≈ besseli(x,1)
    @test legendre_quadrature(dmodified_bessel_1, 7, 0, pi/2) ≈ besseli(x,1)
    @test lobatto_quadrature(dmodified_bessel_1, 8, 0, pi/2) ≈ besseli(x,1)
    @test radau_quadrature(dmodified_bessel_1, 8, 0, pi/2) ≈ besseli(x,1)
    @test rectangle_rule(dmodified_bessel_1, 58999999, 0, pi/2) ≈ besseli(x,1)
    @test simpsons_rule(dmodified_bessel_1, 6, 0, pi/2) ≈ besseli(x,1)
    @test trapezoidal_rule(dmodified_bessel_1, 3, 0, pi/2) ≈ besseli(x,1)
end
