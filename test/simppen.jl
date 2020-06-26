# Simple pendulum test
println("Integrating 1/sqrt(-19.6 sin(x)) from -pi to 0 and comparing the result to the analytical solution of ellipk(1/2)/sqrt(2.45). The integrand has singularities at x = -pi and x=0, so for some of these functions the integration domain has to be itself approximated.")
@testset "Simppen" begin
    println("Running: chebyshev_quadrature with k=1")
    @time @test chebyshev_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 6, 1, -pi, 0) ≈ ellipk(1/2)/sqrt(2.45)
    println("Running: chebyshev_quadrature with k=2")
    @time @test chebyshev_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 7.6983761e7, 2, -pi, 0) ≈ ellipk(1/2)/sqrt(2.45)
    println("Running: chebyshev_quadrature with k=3")
    @time @test chebyshev_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 3.6425654e7, 3, -pi, 0) ≈ ellipk(1/2)/sqrt(2.45)
    println("Running: chebyshev_quadrature with k=4")
    @time @test chebyshev_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 5.0265982e7, 4, -pi, 0) ≈ ellipk(1/2)/sqrt(2.45)
    # Tried up to N=1e7 without Jacobi coming close enough to the analytical solution, and it used a fair chunk of RAM at this N value
    println("Running: jacobi_quadrature with α=β=1. The N value required to get an accurate result is too high to be practical, so only a rough approximation can be arrived at.")
    @time @test abs(jacobi_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 1e6, 1, 1, -pi, 0) - ellipk(1/2)/sqrt(2.45)) < 1e-5
    println("Running: legendre_quadrature.")
    @time @test legendre_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 4.3923464e7, -pi, 0) ≈ ellipk(1/2)/sqrt(2.45)
    println("The following integration functions have to use an approximated domain due to the endpoint singularities.")
    println("Running: lobatto_quadrature on [-pi+1e-6, -1e-6]. Only a rough approximation can be realistically achieved with this function, partly due to the singularities.")
    @time @test abs(lobatto_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 1e6, -pi+1e-6, -1e-6) - ellipk(1/2)/sqrt(2.45)) < 1e-3
    println("Running: radau_quadrature on [-pi+1e-6, -1e-6]. Only a rough approximation can be realistically achieved with this function, partly due to the singularities.")
    @time @test abs(radau_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 1e6, -pi+1e-6, -1e-6) - ellipk(1/2)/sqrt(2.45)) < 1e-3
    println("Running: rectangle_rule on [-pi+1e-8, -1e-8]. Only a rough approximation can be realistically achieved with this function, partly due to the singularities.")
    @time @test abs(rectangle_rule(x -> (-19.6*sin.(x)).^(-0.5), 1e8, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 1e-4
    println("Running: simpsons_rule on [-pi+1e-8, -1e-8]. Only a rough approximation can be realistically achieved with this function, partly due to the singularities.")
    @time @test abs(simpsons_rule(x -> (-19.6*sin.(x)).^(-0.5), 1e8, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 1e-4
    println("Running: trapezoidal_rule on [-pi+1e-8, -1e-8]. Only a rough approximation can be realistically achieved with this function, partly due to the singularities.")
    @time @test abs(trapezoidal_rule(x -> (-19.6*sin.(x)).^(-0.5), 1e8, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 1e-4
end
