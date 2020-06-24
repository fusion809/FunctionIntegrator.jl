# Gaussian curve test
@testset "Gaussian curve" begin
    @test rectangle(x -> exp.(-x.^2), 1e10, 0, 100) ≈ sqrt(pi)/2
    @test simpsons(x -> exp.(-x.^2), 1e4, 0, 100) ≈ sqrt(pi)/2
    @test trapezoidal(x -> exp.(-x.^2), 1e4, 0, 100) ≈ sqrt(pi)/2
    @test hermite_quadrature(x -> 1.0, 1e2, 2) ≈ sqrt(pi)
    @test jacobi_quadrature(x -> exp.(-x.^2), 1e5, 1, 1, 0, 100) ≈ sqrt(pi)/2
    @test laguerre_quadrature(x -> exp.(-x.^2), 1e2, 1) ≈ sqrt(pi)/2
    @test legendre_quadrature(x -> exp.(-x.^2), 1e4, 0, 100) ≈ sqrt(pi)/2
    @test lobatto_quadrature(x -> exp.(-x.^2), 1e4, 0, 100) ≈ sqrt(pi)/2
    @test radau_quadrature(x -> exp.(-x.^2), 1e4, 0, 100) ≈ sqrt(pi)/2
    @test chebyshev_quadrature(x -> exp.(-x.^2), 1e5, 1, 0, 100) ≈ sqrt(pi)/2
    @test chebyshev_quadrature(x -> exp.(-x.^2), 1e5, 2, 0, 100) ≈ sqrt(pi)/2
end
