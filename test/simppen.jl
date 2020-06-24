# Simple pendulum test
@testset "Simple pendulum" begin
    # Singularities exist at the endpoints so rectangle, Simpson's & the trapezoidal rule cannot start and end at them exactly
    @test abs(rectangle(x -> (-19.6*sin.(x)).^(-0.5), 100000000, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 1e-4
    @test abs(simpsons(x -> (-19.6*sin.(x)).^(-0.5), 100000000, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 1e-4
    @test abs(trapezoidal(x -> (-19.6*sin.(x)).^(-0.5), 100000000, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 1e-4
    @test legendre_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 100000000, -pi, 0) ≈ ellipk(1/2)/sqrt(2.45)
    @test chebyshev_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 100, 1, -pi, 0) ≈ ellipk(1/2)/sqrt(2.45)
    # Can't use ≈ for the 2nd Chebyshev quadrature as there's too much error in its estimate for that to work
    @test abs(chebyshev_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 10000, 2, -pi, 0) - ellipk(1/2)/sqrt(2.45)) < 1e-3
end
