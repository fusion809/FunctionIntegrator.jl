# Approximating the integral from -inf to inf of e^(-x^2)/(x^2+1)
sol_9 = exp(1)*pi*erfc(1);

printstyled("Integrating e^(-x^2)/(1+x^2) on the infinite domain [-inf, inf], or failing that on [-100,100], and comparing it to the analytical result exp(1)*pi*erfc(1)\n"; color = :red)
@testset "expnx2datan" begin
    printstyled("Running: adaptive_simpsons_rule with ε=1e-7\n"; color = :magenta)
    @time @test adaptive_simpsons_rule(x -> exp(-x^2)/(x^2+1), -100, 100, 1e-7) ≈ sol_9
    printstyled("Running: chebyshev_quadrature with k=1\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> exp(-x^2)/(x^2+1), 1029, 1, -100, 100) ≈ sol_9
    printstyled("Running: chebyshev_quadrature with k=2\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> exp(-x^2)/(x^2+1), 1028, 2, -100, 100) ≈ sol_9
    printstyled("Running: chebyshev_quadrature with k=3\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> exp(-x^2)/(x^2+1), 514, 3, -100, 100) ≈ sol_9
    printstyled("Running: chebyshev_quadrature with k=4\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> exp(-x^2)/(x^2+1), 514, 4, -100, 100) ≈ sol_9
    printstyled("Running: hermite_quadrature with k=2\n"; color = :magenta)
    @time @test hermite_quadrature(x -> (x^2+1)^(-1), 53, 2) ≈ sol_9
    printstyled("Running: jacobi_quadrature with α=β=1\n"; color = :magenta)
    @time @test jacobi_quadrature(x -> exp(-x^2)/(x^2+1), 1027, 1, 1, -100, 100) ≈ sol_9
    printstyled("Running: laguerre_quadrature with k=1 and multiplying the result by 2 (as Laguerre is only on the semi-infinite domain)\n"; color = :magenta)
    @time @test laguerre_quadrature(x -> 2*exp(-x^2)/(x^2+1), 55, 1) ≈ sol_9
    printstyled("Running: legendre_quadrature\n"; color = :magenta)
    @time @test legendre_quadrature(x -> exp(-x^2)/(x^2+1), 1028, -100, 100) ≈ sol_9
    printstyled("Running: lobatto_quadrature\n"; color = :magenta)
    @time @test lobatto_quadrature(x -> exp(-x^2)/(x^2+1), 1029, -100, 100) ≈ sol_9
    printstyled("Running: radau_quadrature\n"; color = :magenta)
    @time @test radau_quadrature(x -> exp(-x^2)/(x^2+1), 669, -100, 100) ≈ sol_9
    printstyled("Running: rectangle_rule_left\n"; color = :magenta)
    @time @test rectangle_rule_left(x -> exp(-x^2)/(x^2+1), 655, -100, 100) ≈ sol_9
    printstyled("Running: rectangle_rule_midpoint\n"; color = :magenta)
    @time @test rectangle_rule_midpoint(x -> exp(-x^2)/(x^2+1), 655, -100, 100) ≈ sol_9
    printstyled("Running: rectangle_rule_right\n"; color = :magenta)
    @time @test rectangle_rule_right(x -> exp(-x^2)/(x^2+1), 655, -100, 100) ≈ sol_9
    printstyled("Running: rombergs_method\n"; color = :magenta)
    @time @test rombergs_method(x -> exp(-x^2)/(x^2+1), 12, -100, 100) ≈ sol_9
    printstyled("Running: simpsons_rule\n"; color = :magenta)
    @time @test simpsons_rule(x -> exp(-x^2)/(x^2+1), 620, -100, 100) ≈ sol_9
    printstyled("Running: simpsons38_rule\n"; color = :magenta)
    @time @test simpsons38_rule(x -> exp(-x^2)/(x^2+1), 1767, -100, 100) ≈ sol_9
    printstyled("Running: trapezoidal_rule\n"; color = :magenta)
    @time @test trapezoidal_rule(x -> exp(-x^2)/(x^2+1), 655, -100, 100) ≈ sol_9
end
