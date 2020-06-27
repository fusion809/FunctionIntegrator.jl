printstyled("Integrating log(x)/x from 1 to e and comparing the result to the analytical solution of 0.5.\n"; color = :red)
@testset "log(x)/x" begin
    printstyled("Running: chebyshev_quadrature with k=1\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> log.(x).*x.^(-1), 4177, 1, 1, exp(1)) ≈ 0.5
    printstyled("Running: chebyshev_quadrature with k=2\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> log.(x).*x.^(-1), 5906, 2, 1, exp(1)) ≈ 0.5
    printstyled("Running: chebyshev_quadrature with k=3\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> log.(x).*x.^(-1), 4177, 3, 1, exp(1)) ≈ 0.5
    printstyled("Running: chebyshev_quadrature with k=4\n"; color = :magenta)
    @time @test chebyshev_quadrature(x -> log.(x).*x.^(-1), 5907, 4, 1, exp(1)) ≈ 0.5
    printstyled("Running: jacobi_quadrature with α=β=1\n"; color = :magenta)
    @time @test jacobi_quadrature(x -> log.(x).*x.^(-1), 9210, 1, 1, 1, exp(1)) ≈ 0.5
    printstyled("Running: legendre_quadrature\n"; color = :magenta)
    @time @test legendre_quadrature(x -> log.(x).*x.^(-1), 8, 1, exp(1)) ≈ 0.5
    printstyled("Running: lobatto_quadrature\n"; color = :magenta)
    @time @test lobatto_quadrature(x -> log.(x).*x.^(-1), 9, 1, exp(1)) ≈ 0.5
    printstyled("Running: radau_quadrature\n"; color = :magenta)
    @time @test radau_quadrature(x -> log.(x).*x.^(-1), 8, 1, exp(1)) ≈ 0.5
    printstyled("Running: rectangle_rule\n"; color = :magenta)
    @time @test rectangle_rule(x -> log.(x).*x.^(-1), 4.38e7, 1, exp(1)) ≈ 0.5
    printstyled("Running: simpsons_rule\n"; color = :magenta)
    @time @test simpsons_rule(x -> log.(x).*x.^(-1), 100, 1, exp(1)) ≈ 0.5
    printstyled("Running: trapezoidal_rule\n"; color = :magenta)
    @time @test trapezoidal_rule(x -> log.(x).*x.^(-1), 5747, 1, exp(1)) ≈ 0.5
end
