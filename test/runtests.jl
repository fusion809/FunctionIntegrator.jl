using Integration, SpecialFunctions, Test

# Gaussian curve test
@testset "Gaussian curve" begin
    @test rectangle(x -> exp.(-x.^2), 1e10, 0, 100) ≈ sqrt(pi)/2
    @test simpsons(x -> exp.(-x.^2), 10000, 0, 100) ≈ sqrt(pi)/2
    @test trapezoidal(x -> exp.(-x.^2), 10000, 0, 100) ≈ sqrt(pi)/2
    @test legendre_quadrature(x -> exp.(-x.^2), 10000, 0, 100) ≈ sqrt(pi)/2
    @test chebyshev_quadrature(x -> exp.(-x.^2), 1, 100000, 0, 100) ≈ sqrt(pi)/2
    @test chebyshev_quadrature(x -> exp.(-x.^2), 2, 100000, 0, 100) ≈ sqrt(pi)/2
end

# Simple pendulum test
@testset "Simple pendulum" begin
    # Singularities exist at the endpoints so rectangle, Simpson's & the trapezoidal rule cannot start and end at them exactly
    @test abs(rectangle(x -> (-19.6*sin.(x)).^(-0.5), 100000000, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 1e-4
    @test abs(simpsons(x -> (-19.6*sin.(x)).^(-0.5), 100000000, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 1e-4
    @test abs(trapezoidal(x -> (-19.6*sin.(x)).^(-0.5), 100000000, -pi+1e-8, -1e-8) - ellipk(1/2)/sqrt(2.45)) < 1e-4
    @test legendre_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 100000000, -pi, 0) ≈ ellipk(1/2)/sqrt(2.45)
    @test chebyshev_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 1, 100, -pi, 0) ≈ ellipk(1/2)/sqrt(2.45)
    # Can't use ≈ for the 2nd Chebyshev quadrature as there's too much error in its estimate for that to work
    @test abs(chebyshev_quadrature(x -> (-19.6*sin.(x)).^(-0.5), 2, 10000, -pi, 0) - ellipk(1/2)/sqrt(2.45)) < 1e-3
end

# cosine test
@testset "Cosine" begin
    @test abs(rectangle(x -> cos.(x), 1000, 0, pi)) < 1e-14
    @test abs(simpsons(x -> cos.(x), 1000, 0, pi)) < 1e-14
    @test abs(trapezoidal(x -> cos.(x), 1000, 0, pi)) < 1e-14
    @test abs(legendre_quadrature(x -> cos.(x), 1000, 0, pi)) < 1e-14
    @test abs(chebyshev_quadrature(x -> cos.(x), 1, 1000, 0, pi)) < 1e-14
    @test abs(chebyshev_quadrature(x -> cos.(x), 2, 1000, 0, pi)) < 1e-14
end

# logarithm test
@testset "1/x" begin
    @test abs(rectangle(x -> x.^(-1), 1e8, 1, exp(1))-1) < 1e-7
    @test abs(simpsons(x -> x.^(-1), 1000, 1, exp(1))-1) < 1e-10
    @test abs(trapezoidal(x -> x.^(-1), 100000, 1, exp(1))-1) < 1e-10
    @test abs(legendre_quadrature(x -> x.^(-1), 10000, 1, exp(1))-1) < 1e-10
    @test abs(chebyshev_quadrature(x -> x.^(-1), 1, 1000, 1, exp(1))-1) < 1e-5
    @test abs(chebyshev_quadrature(x -> x.^(-1), 2, 1000, 1, exp(1))-1) < 1e-5
end

function cos_cot_fn(x)
    # ^(-1) is used instead of / as / produces a matrix instead of
    # a vector
    (cos.(x).^2).*(cot.(x).+1.0).^(-1)
end

# test an unusual trig integral
@testset "cos_cot_integral" begin
    @test rectangle(cos_cot_fn, 10000, 0, pi/2) ≈ 0.25
    @test simpsons(cos_cot_fn, 10000, 0, pi/2) ≈ 0.25
    @test trapezoidal(cos_cot_fn, 10000, 0, pi/2) ≈ 0.25
    @test legendre_quadrature(cos_cot_fn, 10000, 0, pi/2) ≈ 0.25
    @test chebyshev_quadrature(cos_cot_fn, 1, 10000, 0, pi/2) ≈ 0.25
    @test chebyshev_quadrature(cos_cot_fn, 2, 10000, 0, pi/2) ≈ 0.25
end

@testset "log(x)/x" begin
    @test rectangle(x -> log.(x).*x.^(-1), 1e8, 1, exp(1)) ≈ 0.5
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

@testset "test_7" begin
    @test abs(rectangle(test_7, 1e8, 1, exp(1)) - sol_7) < 1e-8
    @test simpsons(test_7, 10000, 1, exp(1)) ≈ sol_7
    @test trapezoidal(test_7, 100000, 1, exp(1)) ≈ sol_7
    @test legendre_quadrature(test_7, 10000, 1, exp(1)) ≈ sol_7
    @test chebyshev_quadrature(test_7, 1, 10000, 1, exp(1)) ≈ sol_7
    @test chebyshev_quadrature(test_7, 2, 100000, 1, exp(1)) ≈ sol_7
end

# Sinc function
function sinxx(x)
    if x != 0
        return sin.(x).*x.^(-1)
    else
        return 1
    end
end

sol_8 = 1.562225466889056293352345138804502677227824980541083456384;

@testset "sinxx" begin
    @test abs(rectangle(sinxx, 1e8, 0, 100) - sol_8) < 1e-6
    @test simpsons(sinxx, 100000, 0, 100) ≈ sol_8
    @test trapezoidal(sinxx, 100000, 0, 100) ≈ sol_8
    @test legendre_quadrature(sinxx, 10000, 0, 100) ≈ sol_8
    @test chebyshev_quadrature(sinxx, 1, 100000, 0, 100) ≈ sol_8
    @test chebyshev_quadrature(sinxx, 2, 100000, 0, 100) ≈ sol_8
end
