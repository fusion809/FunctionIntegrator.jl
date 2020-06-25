# Test 7
function test_7(x)
    ((x.^3).+1.0).*((x.^4).*(x.+1.0).*((x.^2).+1.0)).^(-1)
end
sol_7 = log(sqrt(2)*exp(1)/(sqrt(exp(2)+1)))+1/2*(exp(-2)-1)+1/3*(1-exp(-3));

@testset "test_7" begin
    @test chebyshev_quadrature(test_7, 8517, 1, 1, exp(1)) ≈ sol_7
    @test chebyshev_quadrature(test_7, 12043, 2, 1, exp(1)) ≈ sol_7
    @test jacobi_quadrature(test_7, 18779, 1, 1, 1, exp(1)) ≈ sol_7
    @test legendre_quadrature(test_7, 10, 1, exp(1)) ≈ sol_7
    @test lobatto_quadrature(test_7, 11, 1, exp(1)) ≈ sol_7
    @test radau_quadrature(test_7, 10, 1, exp(1)) ≈ sol_7
    @test abs(rectangle_rule(test_7, 1e8, 1, exp(1)) - sol_7) < 1e-8
    @test simpsons_rule(test_7, 200, 1, exp(1)) ≈ sol_7
    @test trapezoidal_rule(test_7, 13983, 1, exp(1)) ≈ sol_7
end
