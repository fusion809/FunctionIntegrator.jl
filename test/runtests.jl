using Integration, SpecialFunctions, Test

# Gaussian curve test
N_gauss_test = 10000;
@testset "Gaussian curve" begin
    @test simpsons(x -> exp.(-x.^2), N_gauss_test, 0, 100) ≈ sqrt(pi)/2
    @test trapezoidal(x -> exp.(-x.^2), N_gauss_test, 0, 100) ≈ sqrt(pi)/2
    @test legendre_quadrature(x -> exp.(-x.^2), N_gauss_test, 0, 100) ≈ sqrt(pi)/2
    @test abs(chebyshev_quadrature(x -> exp.(-x.^2), 1, N_gauss_test, 0, 100) - sqrt(pi)/2) < 1e-4
    @test abs(chebyshev_quadrature(x -> exp.(-x.^2), 2, N_gauss_test, 0, 100) - sqrt(pi)/2) < 1e-4
end

# Simple pendulum test
@testset "Simple pendulum" begin
    # Singularities exist at the endpoints so Simpson's & the trapezoidal rule cannot start and end at them exactly
    @test abs(simpsons(x -> (-19.6*sin.(x)).^(-0.5), 100000000, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 1e-4
    @test abs(trapezoidal(x -> (-19.6*sin.(x)).^(-0.5), 100000000, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 1e-4
    @test abs(legendre_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 10000, -pi, 0) - ellipk(1/2)/sqrt(2.45)) < 1e-4
    @test chebyshev_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 1, 100, -pi, 0) ≈ ellipk(1/2)/sqrt(2.45)
    # Can't use ≈ for the 2nd Chebyshev quadrature as there's too much error in its estimate for that to work
    @test abs(chebyshev_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 2, 10000, -pi, 0) - ellipk(1/2)/sqrt(2.45)) < 1e-3
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
@testset "Logarithm" begin
    @test abs(simpsons(x -> x.^(-1), 1000, 1, exp(1))-1) < 1e-10
    @test abs(trapezoidal(x -> x.^(-1), 100000, 1, exp(1))-1) < 1e-10
    @test abs(legendre_quadrature(x -> x.^(-1), 10000, 1, exp(1))-1) < 1e-10
    @test abs(chebyshev_quadrature(x -> x.^(-1), 1, 1000, 1, exp(1))-1) < 1e-5
    @test abs(chebyshev_quadrature(x -> x.^(-1), 2, 1000, 1, exp(1))-1) < 1e-5
end

function cos_cot_fn(x)
    (cos.(x).^2).*(cot.(x).+1.0).^(-1)
end

# test an unusual trig integral
@testset "cos_cot_integral" begin
    @test simpsons(cos_cot_fn, 10000, 0, pi/2) ≈ 0.25
    @test trapezoidal(cos_cot_fn, 10000, 0, pi/2) ≈ 0.25
    @test legendre_quadrature(cos_cot_fn, 10000, 0, pi/2) ≈ 0.25
    @test chebyshev_quadrature(cos_cot_fn, 1, 10000, 0, pi/2) ≈ 0.25
    @test chebyshev_quadrature(cos_cot_fn, 2, 10000, 0, pi/2) ≈ 0.25
end

@testset "log(x)/x" begin
    @test simpsons(x -> log.(x).*x.^(-1), 10000, 1, exp(1)) ≈ 0.5
    @test trapezoidal(x -> log.(x).*x.^(-1), 10000, 1, exp(1)) ≈ 0.5
    @test legendre_quadrature(x -> log.(x).*x.^(-1), 10000, 1, exp(1)) ≈ 0.5
    @test chebyshev_quadrature(x -> log.(x).*x.^(-1), 1, 10000, 1, exp(1)) ≈ 0.5
    @test chebyshev_quadrature(x -> log.(x).*x.^(-1), 2, 10000, 1, exp(1)) ≈ 0.5
end

# Test 7
function test_7(x)
    ((x.^3).+1.0).*((x.^4).*(x.+1.0).*((x.^2).+1.0)).^(-1)
end
sol_7 = log(sqrt(2)*exp(1)/(sqrt(exp(2)+1)))+1/2*(exp(-2)-1)+1/3*(1-exp(-3));

@testset "(x^3+1)/((x^4)(x+1)(x^2+1))" begin
    @test simpsons(test_7, 10000, 1, exp(1)) ≈ sol_7
    @test trapezoidal(test_7, 100000, 1, exp(1)) ≈ sol_7
    @test legendre_quadrature(test_7, 10000, 1, exp(1)) ≈ sol_7
    @test chebyshev_quadrature(test_7, 1, 10000, 1, exp(1)) ≈ sol_7
    @test chebyshev_quadrature(test_7, 2, 100000, 1, exp(1)) ≈ sol_7
end
