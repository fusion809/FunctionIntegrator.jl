# sin(x)/x function
function sinxx(x)
    return sinc(x/pi);
end

sol_8 = 1.562225466889056293352345138804502677227824980541083456384;

printstyled("Integrating sin(x)/x from 0 to 100 and comparing it to the exact result.\n"; color = :red)
@testset "sinxx" begin
    printstyled("Running: chebyshev_quadrature with k=1\n"; color = :magenta)
    @time @test chebyshev_quadrature(sinxx, 29645, 1, 0, 100) ≈ sol_8
    printstyled("Running: chebyshev_quadrature with k=2\n"; color = :magenta)
    @time @test chebyshev_quadrature(sinxx, 41923, 2, 0, 100) ≈ sol_8
    printstyled("Running: chebyshev_quadrature with k=3\n"; color = :magenta)
    @time @test chebyshev_quadrature(sinxx, 42083, 3, 0, 100) ≈ sol_8
    printstyled("Running: chebyshev_quadrature with k=4\n"; color = :magenta)
    @time @test chebyshev_quadrature(sinxx, 29870, 4, 0, 100) ≈ sol_8
    printstyled("Running: jacobi_quadrature with α=β=1\n"; color = :magenta)
    @time @test jacobi_quadrature(sinxx, 65375, 1, 1, 0, 100) ≈ sol_8
    printstyled("Running: legendre_quadrature\n"; color = :magenta)
    @time @test legendre_quadrature(sinxx, 37, 0, 100) ≈ sol_8
    printstyled("Running: lobatto_quadrature\n"; color = :magenta)
    @time @test lobatto_quadrature(sinxx, 38, 0, 100) ≈ sol_8
    printstyled("Running: radau_quadrature\n"; color = :magenta)
    @time @test radau_quadrature(sinxx, 38, 0, 100) ≈ sol_8
    printstyled("Running: rectangle_rule_left; only a rough approximation can be practically achieved using this function\n"; color = :magenta)
    @time @test abs(rectangle_rule_left(sinxx, 1e8, 0, 100) - sol_8) < 5.036e-7
    printstyled("Running: rectangle_rule_midpoint\n"; color = :magenta)
    @time @test rectangle_rule_midpoint(sinxx, 12461, 0, 100) ≈ sol_8
    printstyled("Running: rectangle_rule_right; only a rough approximation can be practically achieved using this function\n"; color = :magenta)
    @time @test abs(rectangle_rule_right(sinxx, 1e8, 0, 100) - sol_8) < 5.015e-7
    printstyled("Running: simpsons_rule\n"; color = :magenta)
    @time @test simpsons_rule(sinxx, 678, 0, 100) ≈ sol_8
    printstyled("Running: trapezoidal_rule\n"; color = :magenta)
    @time @test trapezoidal_rule(sinxx, 17622, 0, 100) ≈ sol_8
end
