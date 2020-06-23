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
function simppen(x)
    (-19.6*sin.(x)).^(-0.5)
end

@testset "Simple pendulum" begin
    @test chebyshev_quadrature(simppen, 1, N_simple_pendulum_test1, -pi, 0) ≈ ellipk(1/2)/sqrt(2.45)
    @test abs(chebyshev_quadrature(simppen, 2, N_simple_pendulum_test2, -pi, 0) - ellipk(1/2)/sqrt(2.45)) < 1e-3
end