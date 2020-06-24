# logarithm test
@testset "1/x" begin
    @test abs(rectangle(x -> x.^(-1), 1e8, 1, exp(1))-1) < 1e-7
    @test abs(simpsons(x -> x.^(-1), 1000, 1, exp(1))-1) < 1e-10
    @test abs(trapezoidal(x -> x.^(-1), 100000, 1, exp(1))-1) < 1e-10
    @test abs(legendre_quadrature(x -> x.^(-1), 10000, 1, exp(1))-1) < 1e-10
    @test abs(chebyshev_quadrature(x -> x.^(-1), 1, 1000, 1, exp(1))-1) < 1e-5
    @test abs(chebyshev_quadrature(x -> x.^(-1), 2, 1000, 1, exp(1))-1) < 1e-5
end
