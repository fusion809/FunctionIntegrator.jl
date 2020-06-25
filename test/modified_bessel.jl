# test an unusual trig integral
x = 1;
function dmodified_bessel_1(theta)
    # ^(-1) is used instead of / as / produces a matrix instead of
    # a vector
    2*x/pi*cosh(x*sin(theta))*(cos(theta))^2
end

@testset "modified_bessel_0" begin
    @test chebyshev_quadrature(dmodified_bessel_1, 4.942e3, 1, 0, pi/2) ≈ besseli(x,1)
    @test chebyshev_quadrature(dmodified_bessel_1, 6.987e3, 2, 0, pi/2) ≈ besseli(x,1)
    @test jacobi_quadrature(dmodified_bessel_1, 1.0896e4, 1, 1, 0, pi/2) ≈ besseli(x,1)
    @test legendre_quadrature(dmodified_bessel_1, 7, 0, pi/2) ≈ besseli(x,1)
    @test lobatto_quadrature(dmodified_bessel_1, 8, 0, pi/2) ≈ besseli(x,1)
    @test radau_quadrature(dmodified_bessel_1, 8, 0, pi/2) ≈ besseli(x,1)
    @test rectangle_rule(dmodified_bessel_1, 5.8999999e7, 0, pi/2) ≈ besseli(x,1)
    @test simpsons_rule(dmodified_bessel_1, 6, 0, pi/2) ≈ besseli(x,1)
    @test trapezoidal_rule(dmodified_bessel_1, 3, 0, pi/2) ≈ besseli(x,1)
end
