# This is the solution that Wolfram Alpha gave for 0 to inf
sol_11 = 0.3592720446932569240312569249695917230704839918143498570462599122401190492471958577982220255972031016;

function sinexpox(x)
    if x!=0
        sin(x^2)*exp(-x)/x
    else
        0
    end
end

printstyled("Integrating sin(x^2)e^(-x)/x from 0 to infinity, with the approximated domain of integration of 0 to 100. The removable singularity at x=0 in the integrand has been removed.\n"; color = :red)
@testset "sinexpox" begin
    printstyled("Running: chebyshev_quadrature with k=1\n"; color = :magenta)
    @time @test chebyshev_quadrature(sinexpox, 655, 1, 0, 100) ≈ sol_11
    printstyled("Running: chebyshev_quadrature with k=2\n"; color = :magenta)
    @time @test chebyshev_quadrature(sinexpox, 654, 2, 0, 100) ≈ sol_11
    printstyled("Running: chebyshev_quadrature with k=3\n"; color = :magenta)
    @time @test chebyshev_quadrature(sinexpox, 666, 3, 0, 100) ≈ sol_11
    printstyled("Running: chebyshev_quadrature with k=4\n"; color = :magenta)
    @time @test chebyshev_quadrature(sinexpox, 576, 4, 0, 100) ≈ sol_11
    printstyled("Running: jacobi_quadrature with α=β=1\n"; color = :magenta)
    @time @test jacobi_quadrature(sinexpox, 598, 1, 1, 0, 100) ≈ sol_11
    printstyled("Running: laguerre_quadrature with k=2\n"; color = :magenta)
    @time @test laguerre_quadrature(x -> sin(x^2)*x^(-1), 4506, 2) ≈ sol_11
    printstyled("Running: legendre_quadrature\n"; color = :magenta)
    @time @test legendre_quadrature(sinexpox, 611, 0, 100) ≈ sol_11
    printstyled("Running: lobatto_quadrature\n"; color = :magenta)
    @time @test lobatto_quadrature(sinexpox, 598, 0, 100) ≈ sol_11
    printstyled("Running: radau_quadrature\n"; color = :magenta)
    @time @test radau_quadrature(sinexpox, 498, 0, 100) ≈ sol_11
    printstyled("Running: rectangle_rule_left\n"; color = :magenta)
    @time @test rectangle_rule_left(sinexpox, 394538, 0, 100) ≈ sol_11
    printstyled("Running: rectangle_rule_midpoint\n"; color = :magenta)
    @time @test rectangle_rule_midpoint(sinexpox, 278981, 0, 100) ≈ sol_11
    printstyled("Running: rectangle_rule_right\n"; color = :magenta)
    @time @test rectangle_rule_right(sinexpox, 394538, 0, 100) ≈ sol_11
    printstyled("Running: simpsons_rule\n"; color = :magenta)
    @time @test simpsons_rule(sinexpox, 4202, 0, 100) ≈ sol_11
    printstyled("Running: simpsons38_rule\n"; color = :magenta)
    @time @test simpsons38_rule(sinexpox, 5148, 0, 100) ≈ sol_11
    printstyled("Running: adaptive_simpsons_rule with ε=1e-8\n"; color = :magenta)
    @time @test adaptive_simpsons_rule(sinexpox, 0, 100, 1e-8) ≈ sol_11
    printstyled("Running: trapezoidal_rule\n"; color = :magenta)
    @time @test trapezoidal_rule(sinexpox, 394538, 0, 100) ≈ sol_11
end
