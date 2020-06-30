# logarithm test
printstyled("Integrating 1/x from 1 to e and comparing the result to the analytical solution, 1.\n"; color = :red)
@testset "1/x" begin
    printstyled("Running: chebyshev_quadrature with k=1\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> x^(-1), 5695, 1, 1, exp(1)) ≈ 1
    printstyled("Running: chebyshev_quadrature with k=2\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> x^(-1), 8053, 2, 1, exp(1)) ≈ 1
    printstyled("Running: chebyshev_quadrature with k=3\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> x^(-1), 6221, 3, 1, exp(1)) ≈ 1
    printstyled("Running: chebyshev_quadrature with k=4\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> x^(-1), 2503, 4, 1, exp(1)) ≈ 1
    printstyled("Running: jacobi_quadrature with α=β=1\n"; color = :magenta)
    @time @test jacobi_quadrature(x -> x^(-1), 12558, 1, 1, 1, exp(1)) ≈ 1
    printstyled("Running: legendre_quadrature\n"; color = :magenta)
    @time @test legendre_quadrature(x -> x^(-1), 7, 1, exp(1)) ≈ 1
    printstyled("Running: lobatto_quadrature\n"; color = :magenta)
    @time @test lobatto_quadrature(x -> x^(-1), 8, 1, exp(1)) ≈ 1
    printstyled("Running: radau_quadrature\n"; color = :magenta)
    @time @test radau_quadrature(x -> x^(-1), 8, 1, exp(1)) ≈ 1
    printstyled("Running: rectangle_rule_left\n"; color = :magenta)
    @time @test rectangle_rule_left(x -> x^(-1), 7.74699e7, 1, exp(1)) ≈ 1
    printstyled("Running: rectangle_rule_midpoint\n"; color = :magenta)
    @time @test rectangle_rule_midpoint(x -> x^(-1), 2672, 1, exp(1)) ≈ 1
    printstyled("Running: rectangle_rule_right\n"; color = :magenta)
    @time @test rectangle_rule_right(x -> x^(-1), 3.6607201e7, 1, exp(1)) ≈ 1
    printstyled("Running: simpsons_rule\n"; color = :magenta)
    @time @test simpsons_rule(x -> x^(-1), 70, 1, exp(1)) ≈ 1
    printstyled("Running: trapezoidal_rule\n"; color = :magenta)
    @time @test trapezoidal_rule(x -> x^(-1), 3779, 1, exp(1)) ≈ 1
end
