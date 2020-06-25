# test an unusual trig integral
function cos_cot_fn(x)
    # ^(-1) is used instead of / as / produces a matrix instead of
    # a vector
    (cos.(x).^2).*(cot.(x).+1.0).^(-1)
end

@testset "coscotint" begin
    @test chebyshev_quadrature(cos_cot_fn, 88, 1, 0, pi/2) ≈ 0.25
    @test chebyshev_quadrature(cos_cot_fn, 90, 2, 0, pi/2) ≈ 0.25
    @test chebyshev_quadrature(cos_cot_fn, 91, 3, 0, pi/2) ≈ 0.25
    @test chebyshev_quadrature(cos_cot_fn, 88, 4, 0, pi/2) ≈ 0.25
    @test jacobi_quadrature(cos_cot_fn, 5, 1, 1, 0, pi/2) ≈ 0.25
    @test legendre_quadrature(cos_cot_fn, 6, 0, pi/2) ≈ 0.25
    @test lobatto_quadrature(cos_cot_fn, 7, 0, pi/2) ≈ 0.25
    @test radau_quadrature(cos_cot_fn, 8, 0, pi/2) ≈ 0.25
    @test rectangle_rule(cos_cot_fn, 7430, 0, pi/2) ≈ 0.25
    @test simpsons_rule(cos_cot_fn, 78, 0, pi/2) ≈ 0.25
    @test trapezoidal_rule(cos_cot_fn, 7430, 0, pi/2) ≈ 0.25
end
