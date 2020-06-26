# Test how accurately the functions integrate an integral
# that should equal I_1 (x)
x = 1;
function dmodified_bessel_1(theta)
    # ^(-1) is used instead of / as / produces a matrix instead of
    # a vector
    2*x/pi*cosh(x*sin(theta))*(cos(theta))^2
end

println("Approximating the modified bessel function I_1(1) and comparing it to the result obtained by SpecialFunctions.")
@testset "modbessel0" begin
    println("Running: chebyshev_quadrature with k=1")
    @time @test chebyshev_quadrature(dmodified_bessel_1, 4942, 1, 0, pi/2) ≈ besseli(x,1)
    println("Running: chebyshev_quadrature with k=2")
    @time @test chebyshev_quadrature(dmodified_bessel_1, 6987, 2, 0, pi/2) ≈ besseli(x,1)
    println("Running: chebyshev_quadrature with k=3")
    @time @test chebyshev_quadrature(dmodified_bessel_1, 6988, 3, 0, pi/2) ≈ besseli(x,1)
    println("Running: chebyshev_quadrature with k=4")
    @time @test chebyshev_quadrature(dmodified_bessel_1, 4941, 4, 0, pi/2) ≈ besseli(x,1)
    println("Running: jacobi_quadrature with α=β=1")
    @time @test jacobi_quadrature(dmodified_bessel_1, 10896, 1, 1, 0, pi/2) ≈ besseli(x,1)
    println("Running: legendre_quadrature")
    @time @test legendre_quadrature(dmodified_bessel_1, 7, 0, pi/2) ≈ besseli(x,1)
    println("Running: lobatto_quadrature")
    @time @test lobatto_quadrature(dmodified_bessel_1, 8, 0, pi/2) ≈ besseli(x,1)
    println("Running: radau_quadrature")
    @time @test radau_quadrature(dmodified_bessel_1, 8, 0, pi/2) ≈ besseli(x,1)
    println("Running: rectangle_rule")
    @time @test rectangle_rule(dmodified_bessel_1, 58999999, 0, pi/2) ≈ besseli(x,1)
    println("Running: simpsons_rule")
    @time @test simpsons_rule(dmodified_bessel_1, 6, 0, pi/2) ≈ besseli(x,1)
    println("Running: trapezoidal_rule")
    @time @test trapezoidal_rule(dmodified_bessel_1, 3, 0, pi/2) ≈ besseli(x,1)
end
