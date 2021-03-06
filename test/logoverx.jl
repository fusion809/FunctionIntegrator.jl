printstyled("Integrating log(x)/x from 1 to e and comparing the result to the analytical solution of 0.5.\n"; color = :red)
@testset "log(x)/x" begin
    printstyled("Running: adaptive_simpsons_rule with ε=1e-7\n"; color = :magenta)
    @time @test adaptive_simpsons_rule(x -> log(x)/x, 1, exp(1), 1e-8) ≈ 0.5
    printstyled("Running: chebyshev_quadrature with k=1\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> log(x)/x, 4177, 1, 1, exp(1)) ≈ 0.5
    printstyled("Running: chebyshev_quadrature with k=2\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> log(x)/x, 5906, 2, 1, exp(1)) ≈ 0.5
    printstyled("Running: chebyshev_quadrature with k=3\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> log(x)/x, 4177, 3, 1, exp(1)) ≈ 0.5
    printstyled("Running: chebyshev_quadrature with k=4\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> log(x)/x, 5907, 4, 1, exp(1)) ≈ 0.5
    printstyled("Running: jacobi_quadrature with α=β=1\n"; color = :magenta)
    @time @test jacobi_quadrature(x -> log(x)/x, 9210, 1, 1, 1, exp(1)) ≈ 0.5
    printstyled("Running: legendre_quadrature\n"; color = :magenta)
    @time @test legendre_quadrature(x -> log(x)/x, 8, 1, exp(1)) ≈ 0.5
    printstyled("Running: lobatto_quadrature\n"; color = :magenta)
    @time @test lobatto_quadrature(x -> log(x)/x, 9, 1, exp(1)) ≈ 0.5
    printstyled("Running: radau_quadrature\n"; color = :magenta)
    @time @test radau_quadrature(x -> log(x)/x, 8, 1, exp(1)) ≈ 0.5
    printstyled("Running: rectangle_rule_left\n"; color = :magenta)
    @time @test rectangle_rule_left(x -> log(x)/x, 4.3800001e7, 1, exp(1)) ≈ 0.5
    printstyled("Running: rectangle_rule_midpoint\n"; color = :magenta)
    @time @test rectangle_rule_midpoint(x -> log(x)/x, 4064, 1, exp(1)) ≈ 0.5
    printstyled("Running: rectangle_rule_right\n"; color = :magenta)
    @time @test rectangle_rule_right(x -> log(x)/x, 4.372011e7, 1, exp(1)) ≈ 0.5
    printstyled("Running: rombergs_method\n"; color = :magenta)
    @time @test rombergs_method(x -> log(x)/x, 5, 1, exp(1)) ≈ 0.5
    printstyled("Running: simpsons_rule\n"; color = :magenta)
    @time @test simpsons_rule(x -> log(x)/x, 50, 1, exp(1)) ≈ 0.5
    printstyled("Running: simpsons38_rule\n"; color = :magenta)
    @time @test simpsons38_rule(x -> log(x)/x, 114, 1, exp(1)) ≈ 0.5
    printstyled("Running: trapezoidal_rule\n"; color = :magenta)
    @time @test trapezoidal_rule(x -> log(x)/x, 5747, 1, exp(1)) ≈ 0.5
end
