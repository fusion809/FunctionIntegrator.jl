# cosine test
println("Integrating cosine from 0 to pi/2 and comparing the result to the analytical result of 1.")
@testset "Cosine" begin
    println("Running: chebyshev_quadrature with k=1")
    @time @test chebyshev_quadrature(x -> cos.(x), 4656, 1, 0, pi/2) ≈ 1
    println("Running: chebyshev_quadrature with k=2")
    @time @test chebyshev_quadrature(x -> cos.(x), 6584, 2, 0, pi/2) ≈ 1
    println("Running: chebyshev_quadrature with k=3")
    @time @test chebyshev_quadrature(x -> cos.(x), 6584, 3, 0, pi/2) ≈ 1
    println("Running: chebyshev_quadrature with k=4")
    @time @test chebyshev_quadrature(x -> cos.(x), 4656, 4, 0, pi/2) ≈ 1
    println("Running: jacobi_quadrature with α=β=1")
    @time @test jacobi_quadrature(x -> cos.(x), 10266, 1, 1, 0, pi/2) ≈ 1
    println("Running: legendre_quadrature")
    @time @test legendre_quadrature(x -> cos.(x), 5, 0, pi/2) ≈ 1
    println("Running: lobatto_quadrature")
    @time @test lobatto_quadrature(x -> cos.(x), 6, 0, pi/2) ≈ 1
    println("Running: radau_quadrature")
    @time @test radau_quadrature(x -> cos.(x), 5, 0, pi/2) ≈ 1
    println("Running: rectangle_rule")
    @time @test rectangle_rule(x -> cos.(x), 5.23705e7, 0, pi/2) ≈ 1
    println("Running: simpsons_rule")
    @time @test simpsons_rule(x -> cos.(x), 40, 0, pi/2) ≈ 1
    println("Running: trapezoidal_rule")
    @time @test trapezoidal_rule(x -> cos.(x), 3715, 0, pi/2) ≈ 1
end
