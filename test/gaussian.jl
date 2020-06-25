# Gaussian curve test
@testset "Gaussian curve" begin
    @test chebyshev_quadrature(x -> exp.(-x.^2), 433, 1, -100, 100) ≈ sqrt(pi)
    @test chebyshev_quadrature(x -> exp.(-x.^2), 432, 2, -100, 100) ≈ sqrt(pi)
    @test hermite_quadrature(x -> 1.0, 1, 2) ≈ sqrt(pi)
    @test jacobi_quadrature(x -> exp.(-x.^2), 87019, 1, 1, -100, 100) ≈ sqrt(pi)
    # Laguerre integrates over [0, inf], so result must be
    # multiplied by 2
    @test laguerre_quadrature(x -> 2*exp.(-x.^2), 46, 1) ≈ sqrt(pi)
    @test legendre_quadrature(x -> exp.(-x.^2), 433, -100, 100) ≈ sqrt(pi)
    @test lobatto_quadrature(x -> exp.(-x.^2), 434, -100, 100) ≈ sqrt(pi)
    @test radau_quadrature(x -> exp.(-x.^2), 349, -100, 100) ≈ sqrt(pi)
    @test rectangle_rule(x -> exp.(-x.^2), 276, -100, 100) ≈ sqrt(pi)
    @test simpsons_rule(x -> exp.(-x.^2), 533, -100, 100) ≈ sqrt(pi)
    @test trapezoidal_rule(x -> exp.(-x.^2), 276, -100, 100) ≈ sqrt(pi)
end
