# Test 7
function partfrac(x)
    (x^3+1.0)/((x^4)*(x+1.0)*(x^2+1.0))
end
sol_7 = log(sqrt(2)*exp(1)/(sqrt(exp(2)+1)))+1/2*(exp(-2)-1)+1/3*(1-exp(-3));

printstyled("Integrating (x^3+1)/(x^4 (x+1)(x^2+1)) from 1 to e and comparing the result to the analytical solution of log(sqrt(2)*exp(1)/(sqrt(exp(2)+1)))+1/2*(exp(-2)-1)+1/3*(1-exp(-3))\n"; color = :red)
@testset "partfrac" begin
    printstyled("Running: adaptive_simpsons_rule with ε=1e-7\n"; color = :magenta)
    @time @test adaptive_simpsons_rule(partfrac, 1, exp(1), 1e-7) ≈ sol_7
    printstyled("Running: chebyshev_quadrature with k=1\n"; color = :magenta)
    @time @test chebyshev_quadrature(partfrac, 8517, 1, 1, exp(1)) ≈ sol_7
    printstyled("Running: chebyshev_quadrature with k=2\n"; color = :magenta)
    @time @test chebyshev_quadrature(partfrac, 12043, 2, 1, exp(1)) ≈ sol_7
    printstyled("Running: chebyshev_quadrature with k=3\n"; color = :magenta)
    @time @test chebyshev_quadrature(partfrac, 11823, 3, 1, exp(1)) ≈ sol_7
    printstyled("Running: chebyshev_quadrature with k=4\n"; color = :magenta)
    @time @test chebyshev_quadrature(partfrac, 8202, 4, 1, exp(1)) ≈ sol_7
    printstyled("Running: jacobi_quadrature with α=β=1\n"; color = :magenta)
    @time @test jacobi_quadrature(partfrac, 18779, 1, 1, 1, exp(1)) ≈ sol_7
    printstyled("Running: legendre_quadrature\n"; color = :magenta)
    @time @test legendre_quadrature(partfrac, 10, 1, exp(1)) ≈ sol_7
    printstyled("Running: lobatto_quadrature\n"; color = :magenta)
    @time @test lobatto_quadrature(partfrac, 11, 1, exp(1)) ≈ sol_7
    printstyled("Running: radau_quadrature\n"; color = :magenta)
    @time @test radau_quadrature(partfrac, 10, 1, exp(1)) ≈ sol_7
    printstyled("Running: rectangle_rule_left; only a rough approximation can be practically achieved using this function.\n"; color = :magenta)
    @time @test abs(rectangle_rule_left(partfrac, 1e8, 1, exp(1)) - sol_7) < 1e-8
    printstyled("Running: rectangle_rule_midpoint\n"; color = :magenta)
    @time @test rectangle_rule_midpoint(partfrac, 9887, 1, exp(1)) ≈ sol_7
    printstyled("Running: rectangle_rule_right; only a rough approximation can be practically achieved using this function.\n"; color = :magenta)
    @time @test abs(rectangle_rule_right(partfrac, 1e8, 1, exp(1)) - sol_7) < 1e-8
    printstyled("Running: simpsons_rule\n"; color = :magenta)
    @time @test simpsons_rule(partfrac, 200, 1, exp(1)) ≈ sol_7
    printstyled("Running: simpsons38_rule\n"; color = :magenta)
    @time @test simpsons38_rule(partfrac, 234, 1, exp(1)) ≈ sol_7
    printstyled("Running: trapezoidal_rule\n"; color = :magenta)
    @time @test trapezoidal_rule(partfrac, 13983, 1, exp(1)) ≈ sol_7
end
