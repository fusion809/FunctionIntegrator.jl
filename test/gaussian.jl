# Gaussian curve test
@testset "Gaussian curve" begin
    @test rectangle(x -> exp.(-x.^2), 1e10, 0, 100) ≈ sqrt(pi)/2
    @test simpsons(x -> exp.(-x.^2), 10000, 0, 100) ≈ sqrt(pi)/2
    @test trapezoidal(x -> exp.(-x.^2), 10000, 0, 100) ≈ sqrt(pi)/2
    @test hermite_quadrature(x -> 1.0, 100, 2) ≈ sqrt(pi)
    @test laguerre_quadrature(x -> exp.(-x.^2), 100, 1) ≈ sqrt(pi)/2
    @test legendre_quadrature(x -> exp.(-x.^2), 10000, 0, 100) ≈ sqrt(pi)/2
    @test chebyshev_quadrature(x -> exp.(-x.^2), 1, 100000, 0, 100) ≈ sqrt(pi)/2
    @test chebyshev_quadrature(x -> exp.(-x.^2), 2, 100000, 0, 100) ≈ sqrt(pi)/2
end
