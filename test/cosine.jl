# cosine test
@testset "Cosine" begin
    @test chebyshev_quadrature(x -> cos.(x), 4656, 1, 0, pi/2) ≈ 1
    @test chebyshev_quadrature(x -> cos.(x), 6584, 2, 0, pi/2) ≈ 1
    @test jacobi_quadrature(x -> cos.(x), 10266, 1, 1, 0, pi/2) ≈ 1
    @test legendre_quadrature(x -> cos.(x), 5, 0, pi/2) ≈ 1
    @test lobatto_quadrature(x -> cos.(x), 6, 0, pi/2) ≈ 1
    @test radau_quadrature(x -> cos.(x), 5, 0, pi/2) ≈ 1
    @test rectangle_rule(x -> cos.(x), 5.23705e7, 0, pi/2) ≈ 1
    @test simpsons_rule(x -> cos.(x), 40, 0, pi/2) ≈ 1
    @test trapezoidal_rule(x -> cos.(x), 3715, 0, pi/2) ≈ 1
end
