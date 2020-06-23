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
simppen_analytical = ellipk(1/2)/sqrt(2.45); 
function simppen(x)
    (-19.6*sin.(x)).^(-0.5)
end

@testset "Simple pendulum" begin
    # Singularities exist at the endpoints so Simpson's & the trapezoidal rule cannot start and end at them exactly
    @test abs(simpsons(simppen, 100000000, -pi+1e-8, -1e-8) - simppen_analytical) < 1e-4
    @test abs(trapezoidal(simppen, 100000000, -pi+1e-8, -1e-8) - simppen_analytical) < 1e-4
    @test abs(legendre_quadrature(simppen, 10000, -pi, 0) - simppen_analytical) < 1e-4
    @test chebyshev_quadrature(simppen, 1, 100, -pi, 0) ≈ simppen_analytical
    # Can't use ≈ for the 2nd Chebyshev quadrature as there's too much error in its estimate for that to work
    @test abs(chebyshev_quadrature(simppen, 2, 10000, -pi, 0) - simppen_analytical) < 1e-3
end

# cosine test
@testset "Cosine" begin
    @test abs(simpsons(x -> cos.(x), 1000, 0, pi)) < 1e-14
    @test abs(trapezoidal(x -> cos.(x), 1000, 0, pi)) < 1e-14
    @test abs(legendre_quadrature(x -> cos.(x), 1000, 0, pi)) < 1e-14
    @test abs(chebyshev_quadrature(x -> cos.(x), 1, 1000, 0, pi)) < 1e-14
    @test abs(chebyshev_quadrature(x -> cos.(x), 2, 1000, 0, pi)) < 1e-14
end

# logarithm test
@testset "Logarithm test" begin
    @test abs(simpsons(x -> x.^(-1), 1000, 1, exp(1))-1) < 1e-10
    @test abs(trapezoidal(x -> x.^(-1), 100000, 1, exp(1))-1) < 1e-10
    @test abs(legendre_quadrature(x -> x.^(-1), 10000, 1, exp(1))-1) < 1e-10
    @test abs(chebyshev_quadrature(x -> x.^(-1), 1, 1000, 1, exp(1))-1) < 1e-5
    @test abs(chebyshev_quadrature(x -> x.^(-1), 2, 1000, 1, exp(1))-1) < 1e-5
end