@testset "log(x)/x" begin
    @test rectangle_rule(x -> log.(x).*x.^(-1), 1e8, 1, exp(1)) ≈ 0.5
    @test simpsons_rule(x -> log.(x).*x.^(-1), 1e4, 1, exp(1)) ≈ 0.5
    @test trapezoidal_rule(x -> log.(x).*x.^(-1), 1e4, 1, exp(1)) ≈ 0.5
    @test jacobi_quadrature(x -> log.(x).*x.^(-1), 1e5, 1, 1, 1, exp(1)) ≈ 0.5
    @test legendre_quadrature(x -> log.(x).*x.^(-1), 1e4, 1, exp(1)) ≈ 0.5
    @test lobatto_quadrature(x -> log.(x).*x.^(-1), 1e4, 1, exp(1)) ≈ 0.5
    @test radau_quadrature(x -> log.(x).*x.^(-1), 1e4, 1, exp(1)) ≈ 0.5
    @test chebyshev_quadrature(x -> log.(x).*x.^(-1), 1e4, 1, 1, exp(1)) ≈ 0.5
    @test chebyshev_quadrature(x -> log.(x).*x.^(-1), 1e4, 2, 1, exp(1)) ≈ 0.5
end
