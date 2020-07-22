# Simple pendulum test
printstyled("Integrating 1/sqrt(-19.6 sin(x)) from -pi to 0 and comparing the result to the analytical solution of ellipk(1/2)/sqrt(2.45) The integrand has singularities at x = -pi and x=0, so for some of these functions the integration domain has to be itself approximated.\n"; color = :red)
@testset "Simppen" begin
    printstyled("Running: adaptive_simpsons_rule on [-pi+1e-8, -1e-8]. Only a rough approximation can be realistically achieved with this function, due to the singularities.\n"; color = :magenta)
    @time @test abs(adaptive_simpsons_rule(x -> (-19.6*sin(x))^(-0.5), -pi+1e-8, -1e-8, 1e-5) - ellipk(1/2)/sqrt(2.45)) < 1e-4
    printstyled("Running: chebyshev_quadrature with k=1\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> (-19.6*sin(x))^(-0.5), 6, 1, -pi, 0) ≈ ellipk(1/2)/sqrt(2.45)
    printstyled("Running: chebyshev_quadrature with k=2\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> (-19.6*sin(x))^(-0.5), 7.6983761e7, 2, -pi, 0) ≈ ellipk(1/2)/sqrt(2.45)
    printstyled("Running: chebyshev_quadrature with k=3\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> (-19.6*sin(x))^(-0.5), 3.6425654e7, 3, -pi, 0) ≈ ellipk(1/2)/sqrt(2.45)
    printstyled("Running: chebyshev_quadrature with k=4\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> (-19.6*sin(x))^(-0.5), 5.0265982e7, 4, -pi, 0) ≈ ellipk(1/2)/sqrt(2.45)
    # Tried up to N=1e7 without Jacobi coming close enough to the analytical solution, and it used a fair chunk of RAM at this N value
    printstyled("Running: jacobi_quadrature with α=β=1. The N value required to get an accurate result is too high to be practical, so only a rough approximation can be arrived at.\n"; color = :magenta)
    @time @test abs(jacobi_quadrature(x -> (-19.6*sin(x))^(-0.5), 1e6, 1, 1, -pi, 0) - ellipk(1/2)/sqrt(2.45)) < 1e-5
    printstyled("Running: legendre_quadrature.\n"; color = :magenta)
    @time @test legendre_quadrature(x -> (-19.6*sin(x))^(-0.5), 4.3923464e7, -pi, 0) ≈ ellipk(1/2)/sqrt(2.45)
    printstyled("The following integration functions have to use an approximated domain due to the endpoint singularities.\n"; color = :magenta)
    printstyled("Running: lobatto_quadrature on [-pi+1e-6, -1e-6]. Only a rough approximation can be realistically achieved with this function, due to the singularities.\n"; color = :magenta)
    @time @test abs(lobatto_quadrature(x -> (-19.6*sin(x))^(-0.5), 1e6, -pi+1e-6, -1e-6) - ellipk(1/2)/sqrt(2.45)) < 1e-3
    printstyled("Running: radau_quadrature on [-pi+1e-6, -1e-6]. Only a rough approximation can be realistically achieved with this function, due to the singularities.\n"; color = :magenta)
    @time @test abs(radau_quadrature(x -> (-19.6*sin(x))^(-0.5), 1e6, -pi+1e-6, -1e-6) - ellipk(1/2)/sqrt(2.45)) < 1e-3
    printstyled("Running: rectangle_rule_left on [-pi+1e-8, -1e-8]. Only a rough approximation can be realistically achieved with this function, due to the singularities.\n"; color = :magenta)
    @time @test abs(rectangle_rule_left(x -> (-19.6*sin(x))^(-0.5), 1e8, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 1e-4
    printstyled("Running: rectangle_rule_midpoint on [-pi+1e-8, -1e-8]. Only a rough approximation can be realistically achieved with this function, due to the singularities.\n"; color = :magenta)
    @time @test abs(rectangle_rule_midpoint(x -> (-19.6*sin(x))^(-0.5), 1e8, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 1.01e-4
    printstyled("Running: rectangle_rule_right on [-pi+1e-8, -1e-8]. Only a rough approximation can be realistically achieved with this function, due to the singularities.\n"; color = :magenta)
    @time @test abs(rectangle_rule_right(x -> (-19.6*sin(x))^(-0.5), 1e8, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 8.649e-5
    printstyled("Running: rombergs_method on [-pi+1e-8, -1e-8]. Only a rough approximation can be realistically achieved with this function, due to the singularities.\n"; color = :magenta)
    @time @test abs(rombergs_method(x -> (-19.6*sin(x))^(-0.5), 24, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 1e-4
    printstyled("Running: simpsons_rule on [-pi+1e-8, -1e-8]. Only a rough approximation can be realistically achieved with this function, due to the singularities.\n"; color = :magenta)
    @time @test abs(simpsons_rule(x -> (-19.6*sin(x))^(-0.5), 1e8, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 1e-4
    printstyled("Running: simpsons38_rule on [-pi+1e-8, -1e-8]. Only a rough approximation can be realistically achieved with this function, due to the singularities.\n"; color = :magenta)
    @time @test abs(simpsons38_rule(x -> (-19.6*sin(x))^(-0.5), 1e8+2, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 1e-4
    printstyled("Running: trapezoidal_rule on [-pi+1e-8, -1e-8]. Only a rough approximation can be realistically achieved with this function, due to the singularities.\n"; color = :magenta)
    @time @test abs(trapezoidal_rule(x -> (-19.6*sin(x))^(-0.5), 1e8, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 1e-4
end