# AiryAi(x)
printstyled("Performing the Airy Ai(x) test, where Ai(x) is integrated on the semi-infinite domain, or an approximation of it, namely [0,100] and the result is compared to the analytical result 1/3.\n"; color = :red)
@testset "Airy Ai(x)" begin
    printstyled("Running: chebyshev_quadrature with k=1\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> airyai(x), 38337, 1, 0, 100) ≈ 1.0/3.0
    printstyled("Running: chebyshev_quadrature with k=2\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> airyai(x), 54215, 2, 0, 100) ≈ 1.0/3.0
    printstyled("Running: chebyshev_quadrature with k=3\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> airyai(x), 54216, 3, 0, 100) ≈ 1.0/3.0
    printstyled("Running: chebyshev_quadrature with k=4\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> airyai(x), 38336, 4, 0, 100) ≈ 1.0/3.0
    printstyled("Running: jacobi_quadrature with α=β=1\n"; color = :magenta)
    @time @test jacobi_quadrature(x -> airyai(x), 84543, 1, 1, 0, 100) ≈ 1.0/3.0
    printstyled("Running: laguerre_quadrature with k=1\n"; color = :magenta)
    @time @test laguerre_quadrature(x -> airyai(x), 15, 1) ≈ 1.0/3.0
    printstyled("Running: legendre_quadrature\n"; color = :magenta)
    @time @test legendre_quadrature(x -> airyai(x), 36, 0, 100) ≈ 1.0/3.0
    printstyled("Running: lobatto_quadrature\n"; color = :magenta)
    @time @test lobatto_quadrature(x -> airyai(x), 37, 0, 100) ≈ 1.0/3.0
    printstyled("Running: radau_quadrature\n"; color = :magenta)
    @time @test radau_quadrature(x -> airyai(x), 37, 0, 100) ≈ 1.0/3.0
    printstyled("Running: rectangle_rule_left. Only a rough approximation can be practically achieved using this function.\n"; color = :magenta)
    @time @test abs(rectangle_rule_left(x -> airyai(x), 1e8, 0, 100) - 1.0/3.0) < 1.775106856505283e-7
    printstyled("Running: rectangle_rule_midpoint.\n"; color = :magenta)
    @time @test rectangle_rule_midpoint(x -> airyai(x), 147348, 0, 100) ≈ 1.0/3.0
    printstyled("Running: rectangle_rule_right. Only a rough approximation can be practically achieved using this function.\n"; color = :magenta)
    @time @test abs(rectangle_rule_right(x -> airyai(x), 1e8, 0, 100) - 1.0/3.0) < 1.775176824944687e-7
    printstyled("Running: simpsons_rule\n"; color = :magenta)
    @time @test simpsons_rule(x -> airyai(x), 2512, 0, 100) ≈ 1.0/3.0
    printstyled("Running: simpsons38_rule\n"; color = :magenta)
    @time @test simpsons38_rule(x -> airyai(x), 3075, 0, 100) ≈ 1.0/3.0
    printstyled("Running: adaptive_simpsons_rule\n"; color = :magenta)
    @time @test adaptive_simpsons_rule(x -> airyai(x), 0, 100, 1e-8) ≈ 1.0/3.0
    printstyled("Running: trapezoidal_rule\n"; color = :magenta)
    @time @test trapezoidal_rule(x -> airyai(x), 208381, 0, 100) ≈ 1.0/3.0
end
