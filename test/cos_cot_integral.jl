# test an unusual trig integral
function cos_cot_fn(x)
    # ^(-1) is used instead of / as / produces a matrix instead of
    # a vector
    (cos.(x).^2).*(cot.(x).+1.0).^(-1)
end

printstyled("Integrating cos^2(x)/(1+cot(x)) from 0 to pi/2 and comparing the results to the analytical result 0.25.\n"; color = :red)
@testset "coscotint" begin
    printstyled("Running: chebyshev_quadrature with k=1\n"; color = :magenta)
    @time @test chebyshev_quadrature(cos_cot_fn, 88, 1, 0, pi/2) ≈ 0.25
    printstyled("Running: chebyshev_quadrature with k=2\n"; color = :magenta)
    @time @test chebyshev_quadrature(cos_cot_fn, 90, 2, 0, pi/2) ≈ 0.25
    printstyled("Running: chebyshev_quadrature with k=3\n"; color = :magenta)
    @time @test chebyshev_quadrature(cos_cot_fn, 91, 3, 0, pi/2) ≈ 0.25
    printstyled("Running: chebyshev_quadrature with k=4\n"; color = :magenta)
    @time @test chebyshev_quadrature(cos_cot_fn, 88, 4, 0, pi/2) ≈ 0.25
    printstyled("Running: jacobi_quadrature with α=β=1\n"; color = :magenta)
    @time @test jacobi_quadrature(cos_cot_fn, 5, 1, 1, 0, pi/2) ≈ 0.25
    printstyled("Running: legendre_quadrature\n"; color = :magenta)
    @time @test legendre_quadrature(cos_cot_fn, 6, 0, pi/2) ≈ 0.25
    printstyled("Running: lobatto_quadrature\n"; color = :magenta)
    @time @test lobatto_quadrature(cos_cot_fn, 7, 0, pi/2) ≈ 0.25
    printstyled("Running: radau_quadrature\n"; color = :magenta)
    @time @test radau_quadrature(cos_cot_fn, 8, 0, pi/2) ≈ 0.25
    printstyled("Running: rectangle_rule_left\n"; color = :magenta)
    @time @test rectangle_rule_left(cos_cot_fn, 7430, 0, pi/2) ≈ 0.25
    printstyled("Running: rectangle_rule_midpoint\n"; color = :magenta)
    @time @test rectangle_rule_midpoint(cos_cot_fn, 5254, 0, pi/2) ≈ 0.25
    printstyled("Running: rectangle_rule_right\n"; color = :magenta)
    @time @test rectangle_rule_right(cos_cot_fn, 7430, 0, pi/2) ≈ 0.25
    printstyled("Running: simpsons_rule\n"; color = :magenta)
    @time @test simpsons_rule(cos_cot_fn, 78, 0, pi/2) ≈ 0.25
    printstyled("Running: trapezoidal_rule\n"; color = :magenta)
    @time @test trapezoidal_rule(cos_cot_fn, 7430, 0, pi/2) ≈ 0.25
end
