println("Integrating log(x)/x from 1 to e and comparing the result to the analytical solution of 0.5.")
@testset "log(x)/x" begin
    println("Running: chebyshev_quadrature with k=1")
    @time @test chebyshev_quadrature(x -> log.(x).*x.^(-1), 4177, 1, 1, exp(1)) ≈ 0.5
    println("Running: chebyshev_quadrature with k=2")
    @time @test chebyshev_quadrature(x -> log.(x).*x.^(-1), 5906, 2, 1, exp(1)) ≈ 0.5
    println("Running: chebyshev_quadrature with k=3")
    @time @test chebyshev_quadrature(x -> log.(x).*x.^(-1), 4177, 3, 1, exp(1)) ≈ 0.5
    println("Running: chebyshev_quadrature with k=4")
    @time @test chebyshev_quadrature(x -> log.(x).*x.^(-1), 5907, 4, 1, exp(1)) ≈ 0.5
    println("Running: jacobi_quadrature with α=β=1")
    @time @test jacobi_quadrature(x -> log.(x).*x.^(-1), 9210, 1, 1, 1, exp(1)) ≈ 0.5
    println("Running: legendre_quadrature")
    @time @test legendre_quadrature(x -> log.(x).*x.^(-1), 8, 1, exp(1)) ≈ 0.5
    println("Running: lobatto_quadrature")
    @time @test lobatto_quadrature(x -> log.(x).*x.^(-1), 9, 1, exp(1)) ≈ 0.5
    println("Running: radau_quadrature")
    @time @test radau_quadrature(x -> log.(x).*x.^(-1), 8, 1, exp(1)) ≈ 0.5
    println("Running: rectangle_rule")
    @time @test rectangle_rule(x -> log.(x).*x.^(-1), 4.38e7, 1, exp(1)) ≈ 0.5
    println("Running: simpsons_rule")
    @time @test simpsons_rule(x -> log.(x).*x.^(-1), 100, 1, exp(1)) ≈ 0.5
    println("Running: trapezoidal_rule")
    @time @test trapezoidal_rule(x -> log.(x).*x.^(-1), 5747, 1, exp(1)) ≈ 0.5
end
