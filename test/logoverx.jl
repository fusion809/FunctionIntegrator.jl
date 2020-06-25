@testset "log(x)/x" begin
    @test chebyshev_quadrature(x -> log.(x).*x.^(-1), 4177, 1, 1, exp(1)) ≈ 0.5
    @test chebyshev_quadrature(x -> log.(x).*x.^(-1), 5906, 2, 1, exp(1)) ≈ 0.5
    @test chebyshev_quadrature(x -> log.(x).*x.^(-1), 4177, 3, 1, exp(1)) ≈ 0.5
    @test chebyshev_quadrature(x -> log.(x).*x.^(-1), 5907, 4, 1, exp(1)) ≈ 0.5
    @test jacobi_quadrature(x -> log.(x).*x.^(-1), 9210, 1, 1, 1, exp(1)) ≈ 0.5
    @test legendre_quadrature(x -> log.(x).*x.^(-1), 8, 1, exp(1)) ≈ 0.5
    @test lobatto_quadrature(x -> log.(x).*x.^(-1), 9, 1, exp(1)) ≈ 0.5
    @test radau_quadrature(x -> log.(x).*x.^(-1), 8, 1, exp(1)) ≈ 0.5
    @test rectangle_rule(x -> log.(x).*x.^(-1), 4.38e7, 1, exp(1)) ≈ 0.5
    @test simpsons_rule(x -> log.(x).*x.^(-1), 100, 1, exp(1)) ≈ 0.5
    @test trapezoidal_rule(x -> log.(x).*x.^(-1), 5747, 1, exp(1)) ≈ 0.5
end
