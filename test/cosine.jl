# cosine test
printstyled("Integrating cosine from 0 to pi/2 and comparing the result to the analytical result of 1.\n"; color = :red)
@testset "Cosine" begin
    printstyled("Running: chebyshev_quadrature with k=1\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> cos(x), 4656, 1, 0, pi/2) ≈ 1
    printstyled("Running: chebyshev_quadrature with k=2\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> cos(x), 6584, 2, 0, pi/2) ≈ 1
    printstyled("Running: chebyshev_quadrature with k=3\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> cos(x), 6584, 3, 0, pi/2) ≈ 1
    printstyled("Running: chebyshev_quadrature with k=4\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> cos(x), 4656, 4, 0, pi/2) ≈ 1
    printstyled("Running: jacobi_quadrature with α=β=1\n"; color = :magenta)
    @time @test jacobi_quadrature(x -> cos(x), 10266, 1, 1, 0, pi/2) ≈ 1
    printstyled("Running: legendre_quadrature\n"; color = :magenta)
    @time @test legendre_quadrature(x -> cos(x), 5, 0, pi/2) ≈ 1
    printstyled("Running: lobatto_quadrature\n"; color = :magenta)
    @time @test lobatto_quadrature(x -> cos(x), 6, 0, pi/2) ≈ 1
    printstyled("Running: radau_quadrature\n"; color = :magenta)
    @time @test radau_quadrature(x -> cos(x), 5, 0, pi/2) ≈ 1
    printstyled("Running: rectangle_rule_left\n"; color = :magenta)
    @time @test rectangle_rule_left(x -> cos(x), 5.23705e7, 0, pi/2) ≈ 1
    printstyled("Running: rectangle_rule_midpoint\n"; color = :magenta)
    @time @test rectangle_rule_midpoint(x -> cos(x), 2627, 0, pi/2) ≈ 1
    printstyled("Running: rectangle_rule_right\n"; color = :magenta)
    @time @test rectangle_rule_right(x -> cos(x), 5.2441522e7, 0, pi/2) ≈ 1
    printstyled("Running: simpsons_rule\n"; color = :magenta)
    @time @test simpsons_rule(x -> cos(x), 40, 0, pi/2) ≈ 1
    printstyled("Running: trapezoidal_rule\n"; color = :magenta)
    @time @test trapezoidal_rule(x -> cos(x), 3715, 0, pi/2) ≈ 1
end
