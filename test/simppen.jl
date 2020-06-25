# Simple pendulum test
@testset "Simple pendulum" begin
    # Singularities exist at the endpoints so rectangle_rule, Simpson's & the trapezoidal_rule rule cannot start and end at them exactly
    @test abs(rectangle_rule(x -> (-19.6*sin.(x)).^(-0.5), 1e8, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 1e-4
    @test abs(simpsons_rule(x -> (-19.6*sin.(x)).^(-0.5), 1e8, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 1e-4
    @test abs(trapezoidal_rule(x -> (-19.6*sin.(x)).^(-0.5), 1e8, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 1e-4
    @test abs(jacobi_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 1e6, 1, 1, -pi, 0) - ellipk(1/2)/sqrt(2.45)) < 1e-5
    @test legendre_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 1e8, -pi, 0) ≈ ellipk(1/2)/sqrt(2.45)
    @test abs(lobatto_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 1e6, -pi+1e-6, -1e-6) - ellipk(1/2)/sqrt(2.45)) < 1e-3
    @test abs(radau_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 1e6, -pi+1e-6, -1e-6) - ellipk(1/2)/sqrt(2.45)) < 1e-3
    @test chebyshev_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 1e2, 1, -pi, 0) ≈ ellipk(1/2)/sqrt(2.45)
    # Can't use ≈ for the 2nd Chebyshev quadrature as there's too much error in its estimate for that to work
    @test abs(chebyshev_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 1e4, 2, -pi, 0) - ellipk(1/2)/sqrt(2.45)) < 1e-3
end
