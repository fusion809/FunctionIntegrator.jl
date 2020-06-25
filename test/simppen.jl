# Simple pendulum test
@testset "Simppen" begin
    @test chebyshev_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 6, 1, -pi, 0) ≈ ellipk(1/2)/sqrt(2.45)
    @test chebyshev_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 7.6983761e7, 2, -pi, 0) ≈ ellipk(1/2)/sqrt(2.45)
    @test chebyshev_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 3.6425654e7, 3, -pi, 0) ≈ ellipk(1/2)/sqrt(2.45)
    @test chebyshev_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 5.0265982e7, 4, -pi, 0) ≈ ellipk(1/2)/sqrt(2.45)
    # Tried up to N=1e7 without Jacobi coming close enough to the analytical solution, and it used a fair chunk of RAM at this N value
    @test abs(jacobi_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 1e6, 1, 1, -pi, 0) - ellipk(1/2)/sqrt(2.45)) < 1e-5
    @test legendre_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 4.3923464e7, -pi, 0) ≈ ellipk(1/2)/sqrt(2.45)
    @test abs(lobatto_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 1e6, -pi+1e-6, -1e-6) - ellipk(1/2)/sqrt(2.45)) < 1e-3
    @test abs(radau_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 1e6, -pi+1e-6, -1e-6) - ellipk(1/2)/sqrt(2.45)) < 1e-3
    # Singularities exist at the endpoints so rectangle_rule, Simpson's & the trapezoidal_rule rule cannot start and end at them exactly
    @test abs(rectangle_rule(x -> (-19.6*sin.(x)).^(-0.5), 1e8, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 1e-4
    @test abs(simpsons_rule(x -> (-19.6*sin.(x)).^(-0.5), 1e8, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 1e-4
    @test abs(trapezoidal_rule(x -> (-19.6*sin.(x)).^(-0.5), 1e8, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 1e-4
end
