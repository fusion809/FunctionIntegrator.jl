# Test 7
function test_7(x)
    ((x.^3).+1.0).*((x.^4).*(x.+1.0).*((x.^2).+1.0)).^(-1)
end
sol_7 = log(sqrt(2)*exp(1)/(sqrt(exp(2)+1)))+1/2*(exp(-2)-1)+1/3*(1-exp(-3));

@testset "test_7" begin
    @test abs(rectangle(test_7, 1e8, 1, exp(1)) - sol_7) < 1e-8
    @test simpsons(test_7, 10000, 1, exp(1)) ≈ sol_7
    @test trapezoidal(test_7, 100000, 1, exp(1)) ≈ sol_7
    @test legendre_quadrature(test_7, 10000, 1, exp(1)) ≈ sol_7
    @test chebyshev_quadrature(test_7, 1, 10000, 1, exp(1)) ≈ sol_7
    @test chebyshev_quadrature(test_7, 2, 100000, 1, exp(1)) ≈ sol_7
end
