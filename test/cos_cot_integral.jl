# test an unusual trig integral
function cos_cot_fn(x)
    # ^(-1) is used instead of / as / produces a matrix instead of
    # a vector
    (cos.(x).^2).*(cot.(x).+1.0).^(-1)
end

@testset "cos_cot_integral" begin
    @test rectangle(cos_cot_fn, 1e4, 0, pi/2) ≈ 0.25
    @test simpsons(cos_cot_fn, 1e4, 0, pi/2) ≈ 0.25
    @test trapezoidal(cos_cot_fn, 1e4, 0, pi/2) ≈ 0.25
    @test jacobi_quadrature(cos_cot_fn, 1e4, 1, 1, 0, pi/2) ≈ 0.25
    @test legendre_quadrature(cos_cot_fn, 1e4, 0, pi/2) ≈ 0.25
    @test lobatto_quadrature(cos_cot_fn, 1e4, 0, pi/2) ≈ 0.25
    @test chebyshev_quadrature(cos_cot_fn, 1e4, 1, 0, pi/2) ≈ 0.25
    @test chebyshev_quadrature(cos_cot_fn, 1e4, 2, 0, pi/2) ≈ 0.25
    @test radau_quadrature(cos_cot_fn, 1e4, 0, pi/2) ≈ 0.25
end
