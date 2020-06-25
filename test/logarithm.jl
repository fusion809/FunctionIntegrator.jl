# logarithm test
@testset "1/x" begin
    @test abs(chebyshev_quadrature(x -> x.^(-1), 1e3, 1, 1, exp(1))-1) < 1e-5
    @test abs(chebyshev_quadrature(x -> x.^(-1), 1e3, 2, 1, exp(1))-1) < 1e-5
    @test abs(jacobi_quadrature(x -> x.^(-1), 1e5, 1, 1, 1, exp(1))-1) < 1e-9
    @test abs(legendre_quadrature(x -> x.^(-1), 1e4, 1, exp(1))-1) < 1e-10
    @test abs(lobatto_quadrature(x -> x.^(-1), 1e4, 1, exp(1))-1) < 1e-10
    @test abs(radau_quadrature(x -> x.^(-1), 1e4, 1, exp(1))-1) < 1e-10
    @test abs(rectangle_rule(x -> x.^(-1), 1e8, 1, exp(1))-1) < 1e-7
    @test abs(simpsons_rule(x -> x.^(-1), 1e3, 1, exp(1))-1) < 1e-10
    @test abs(trapezoidal_rule(x -> x.^(-1), 1e5, 1, exp(1))-1) < 1e-10
end
