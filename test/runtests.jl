using Integration, SpecialFunctions, Test

# Gaussian curve test
N_gauss_test = 10000;
function expnx2(x)
    exp.(-x.^2)
end

@testset "Gaussian curve" begin
    @test simpsons(expnx2, N_gauss_test, 0, 100) ≈ sqrt(pi)/2
    @test trapezoidal(expnx2, N_gauss_test, 0, 100) ≈ sqrt(pi)/2
    @test legendre_quadrature(expnx2, N_gauss_test, 0, 100) ≈ sqrt(pi)/2
    @test abs(chebyshev_quadrature(expnx2, 1, N_gauss_test, 0, 100) - sqrt(pi)/2) < 1e-4
    @test abs(chebyshev_quadrature(expnx2, 2, N_gauss_test, 0, 100) - sqrt(pi)/2) < 1e-4
end

# Simple pendulum test
N_simple_pendulum_test1 = 100;
N_simple_pendulum_test2 = 10000;
simppen_analytical = ellipk(1/2)/sqrt(2.45); 
function simppen(x)
    (-19.6*sin.(x)).^(-0.5)
end

@testset "Simple pendulum" begin
    # Singularities exist at the endpoints so Simpson's & the trapezoidal rule cannot start and end at them exactly
    @test abs(simpsons(simppen, 100000000, -pi+1e-8, -1e-8) - simppen_analytical) < 1e-4
    @test abs(trapezoidal(simppen, 100000000, -pi+1e-8, -1e-8) - simppen_analytical) < 1e-4
    @test abs(legendre_quadrature(simppen, 100000000, -pi, 0) - simppen_analytical) < 1e-4
    @test chebyshev_quadrature(simppen, 1, N_simple_pendulum_test1, -pi, 0) ≈ simppen_analytical
    # Can't use ≈ for the 2nd Chebyshev quadrature as there's too much error in its estimate for that to work
    @test abs(chebyshev_quadrature(simppen, 2, N_simple_pendulum_test2, -pi, 0) - simppen_analytical) < 1e-3
end