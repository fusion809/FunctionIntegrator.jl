# Approximating the integral from -inf to inf of e^(-x^2)/(x^2+1)
sol_9 = exp(1)*pi*erfc(1);

printstyled("Integrating e^(-x^2)/(1+x^2) on the infinite domain [-inf, inf], or failing that on [-100,100], and comparing it to the analytical result exp(1)*pi*erfc(1).\n"; color = :red)
@testset "expnx2datan" begin
    printstyled("Running: chebyshev_quadrature with k=1\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> exp.(-x.^2).*(x.^2+1).^(-1), 1029, 1, -100, 100) ≈ sol_9
    printstyled("Running: chebyshev_quadrature with k=2\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> exp.(-x.^2).*(x.^2+1).^(-1), 1028, 2, -100, 100) ≈ sol_9
    printstyled("Running: chebyshev_quadrature with k=3\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> exp.(-x.^2).*(x.^2+1).^(-1), 514, 3, -100, 100) ≈ sol_9
    printstyled("Running: chebyshev_quadrature with k=4\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> exp.(-x.^2).*(x.^2+1).^(-1), 514, 4, -100, 100) ≈ sol_9
    @time @test hermite_quadrature(x -> (x.^2+1).^(-1), 53, 2) ≈ sol_9
    printstyled("Running: jacobi_quadrature with α=β=1\n"; color = :magenta)
    @time @test jacobi_quadrature(x -> exp.(-x.^2).*(x.^2+1).^(-1), 1027, 1, 1, -100, 100) ≈ sol_9
    @time @test laguerre_quadrature(x -> 2*exp.(-x.^2).*(x.^2+1).^(-1), 55, 1) ≈ sol_9
    printstyled("Running: legendre_quadrature\n"; color = :magenta)
    @time @test legendre_quadrature(x -> exp.(-x.^2).*(x.^2+1).^(-1), 1028, -100, 100) ≈ sol_9
    printstyled("Running: lobatto_quadrature\n"; color = :magenta)
    @time @test lobatto_quadrature(x -> exp.(-x.^2).*(x.^2+1).^(-1), 1029, -100, 100) ≈ sol_9
    printstyled("Running: radau_quadrature\n"; color = :magenta)
    @time @test radau_quadrature(x -> exp.(-x.^2).*(x.^2+1).^(-1), 669, -100, 100) ≈ sol_9
    printstyled("Running: rectangle_rule\n"; color = :magenta)
    @time @test rectangle_rule(x -> exp.(-x.^2).*(x.^2+1).^(-1), 655, -100, 100) ≈ sol_9
    printstyled("Running: simpsons_rule\n"; color = :magenta)
    @time @test simpsons_rule(x -> exp.(-x.^2).*(x.^2+1).^(-1), 1233, -100, 100) ≈ sol_9
    printstyled("Running: trapezoidal_rule\n"; color = :magenta)
    @time @test trapezoidal_rule(x -> exp.(-x.^2).*(x.^2+1).^(-1), 655, -100, 100) ≈ sol_9
end
