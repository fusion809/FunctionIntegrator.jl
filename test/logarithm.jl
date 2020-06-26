# logarithm test
println("Integrating 1/x from 1 to e and comparing the result to the analytical solution, 1.")
@testset "1/x" begin
    println("Running: chebyshev_quadrature with k=1")
    @time @test chebyshev_quadrature(x -> x.^(-1), 5695, 1, 1, exp(1)) ≈ 1
    println("Running: chebyshev_quadrature with k=2")
    @time @test chebyshev_quadrature(x -> x.^(-1), 8053, 2, 1, exp(1)) ≈ 1
    println("Running: chebyshev_quadrature with k=3")
    @time @test chebyshev_quadrature(x -> x.^(-1), 6221, 3, 1, exp(1)) ≈ 1
    println("Running: chebyshev_quadrature with k=4")
    @time @test chebyshev_quadrature(x -> x.^(-1), 2503, 4, 1, exp(1)) ≈ 1
    println("Running: jacobi_quadrature with α=β=1")
    @time @test jacobi_quadrature(x -> x.^(-1), 12558, 1, 1, 1, exp(1)) ≈ 1
    println("Running: legendre_quadrature")
    @time @test legendre_quadrature(x -> x.^(-1), 7, 1, exp(1)) ≈ 1
    println("Running: lobatto_quadrature")
    @time @test lobatto_quadrature(x -> x.^(-1), 8, 1, exp(1)) ≈ 1
    println("Running: radau_quadrature")
    @time @test radau_quadrature(x -> x.^(-1), 8, 1, exp(1)) ≈ 1
    println("Running: rectangle_rule")
    @time @test rectangle_rule(x -> x.^(-1), 7.74699e7, 1, exp(1)) ≈ 1
    println("Running: simpsons_rule")
    @time @test simpsons_rule(x -> x.^(-1), 70, 1, exp(1)) ≈ 1
    println("Running: trapezoidal_rule")
    @time @test trapezoidal_rule(x -> x.^(-1), 3779, 1, exp(1)) ≈ 1
end
