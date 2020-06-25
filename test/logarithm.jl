# logarithm test
@testset "1/x" begin
    @test chebyshev_quadrature(x -> x.^(-1), 5695, 1, 1, exp(1)) ≈ 1
    @test chebyshev_quadrature(x -> x.^(-1), 8053, 2, 1, exp(1)) ≈ 1
    @test chebyshev_quadrature(x -> x.^(-1), 6221, 3, 1, exp(1)) ≈ 1
    @test chebyshev_quadrature(x -> x.^(-1), 2503, 4, 1, exp(1)) ≈ 1
    @test jacobi_quadrature(x -> x.^(-1), 12558, 1, 1, 1, exp(1)) ≈ 1
    @test legendre_quadrature(x -> x.^(-1), 7, 1, exp(1)) ≈ 1
    @test lobatto_quadrature(x -> x.^(-1), 8, 1, exp(1)) ≈ 1
    @test radau_quadrature(x -> x.^(-1), 8, 1, exp(1)) ≈ 1
    @test rectangle_rule(x -> x.^(-1), 7.74699e7, 1, exp(1)) ≈ 1
    @test simpsons_rule(x -> x.^(-1), 70, 1, exp(1)) ≈ 1
    @test trapezoidal_rule(x -> x.^(-1), 3779, 1, exp(1)) ≈ 1
end
