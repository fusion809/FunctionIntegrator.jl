# test an unusual trig integral
function cos_cot_fn(x)
    # ^(-1) is used instead of / as / produces a matrix instead of
    # a vector
    (cos.(x).^2).*(cot.(x).+1.0).^(-1)
end

println("Integrating cos^2(x)/(1+cot(x)) from 0 to pi/2 and comparing the results to the analytical result 0.25.")
@testset "coscotint" begin
    println("Running: chebyshev_quadrature with k=1")
    @time @test chebyshev_quadrature(cos_cot_fn, 88, 1, 0, pi/2) ≈ 0.25
    println("Running: chebyshev_quadrature with k=2")
    @time @test chebyshev_quadrature(cos_cot_fn, 90, 2, 0, pi/2) ≈ 0.25
    println("Running: chebyshev_quadrature with k=3")
    @time @test chebyshev_quadrature(cos_cot_fn, 91, 3, 0, pi/2) ≈ 0.25
    println("Running: chebyshev_quadrature with k=4")
    @time @test chebyshev_quadrature(cos_cot_fn, 88, 4, 0, pi/2) ≈ 0.25
    println("Running: jacobi_quadrature with α=β=1")
    @time @test jacobi_quadrature(cos_cot_fn, 5, 1, 1, 0, pi/2) ≈ 0.25
    println("Running: legendre_quadrature")
    @time @test legendre_quadrature(cos_cot_fn, 6, 0, pi/2) ≈ 0.25
    println("Running: lobatto_quadrature")
    @time @test lobatto_quadrature(cos_cot_fn, 7, 0, pi/2) ≈ 0.25
    println("Running: radau_quadrature")
    @time @test radau_quadrature(cos_cot_fn, 8, 0, pi/2) ≈ 0.25
    println("Running: rectangle_rule")
    @time @test rectangle_rule(cos_cot_fn, 7430, 0, pi/2) ≈ 0.25
    println("Running: simpsons_rule")
    @time @test simpsons_rule(cos_cot_fn, 78, 0, pi/2) ≈ 0.25
    println("Running: trapezoidal_rule")
    @time @test trapezoidal_rule(cos_cot_fn, 7430, 0, pi/2) ≈ 0.25
end
