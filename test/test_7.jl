# Test 7
function test_7(x)
    ((x.^3).+1.0).*((x.^4).*(x.+1.0).*((x.^2).+1.0)).^(-1)
end
sol_7 = log(sqrt(2)*exp(1)/(sqrt(exp(2)+1)))+1/2*(exp(-2)-1)+1/3*(1-exp(-3));

@testset "test_7" begin
    println("Running: chebyshev_quadrature with k=1")
    @time @test chebyshev_quadrature(test_7, 8517, 1, 1, exp(1)) ≈ sol_7
    println("Running: chebyshev_quadrature with k=2")
    @time @test chebyshev_quadrature(test_7, 12043, 2, 1, exp(1)) ≈ sol_7
    println("Running: chebyshev_quadrature with k=3")
    @time @test chebyshev_quadrature(test_7, 11823, 3, 1, exp(1)) ≈ sol_7
    println("Running: chebyshev_quadrature with k=4")
    @time @test chebyshev_quadrature(test_7, 8202, 4, 1, exp(1)) ≈ sol_7
    println("Running: jacobi_quadrature with α=β=1")
    @time @test jacobi_quadrature(test_7, 18779, 1, 1, 1, exp(1)) ≈ sol_7
    println("Running: legendre_quadrature")
    @time @test legendre_quadrature(test_7, 10, 1, exp(1)) ≈ sol_7
    println("Running: lobatto_quadrature")
    @time @test lobatto_quadrature(test_7, 11, 1, exp(1)) ≈ sol_7
    println("Running: radau_quadrature")
    @time @test radau_quadrature(test_7, 10, 1, exp(1)) ≈ sol_7
    println("Running: rectangle_rule; only a rough approximation can be practically achieved using this function.")
    @time @test abs(rectangle_rule(test_7, 1e8, 1, exp(1)) - sol_7) < 1e-8
    println("Running: simpsons_rule")
    @time @test simpsons_rule(test_7, 200, 1, exp(1)) ≈ sol_7
    println("Running: trapezoidal_rule")
    @time @test trapezoidal_rule(test_7, 13983, 1, exp(1)) ≈ sol_7
end
