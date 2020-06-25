# cosine test
@testset "Cosine" begin
    @test abs(rectangle_rule(x -> cos.(x), 1e3, 0, pi)) < 1e-14
    @test abs(simpsons_rule(x -> cos.(x), 1e3, 0, pi)) < 1e-14
    @test abs(trapezoidal_rule(x -> cos.(x), 1e3, 0, pi)) < 1e-14
    @test abs(jacobi_quadrature(x -> cos.(x), 1e4, 1, 1, 0, pi)) < 1e-14
    @test abs(legendre_quadrature(x -> cos.(x), 1e4, 0, pi)) < 1e-14
    @test abs(lobatto_quadrature(x -> cos.(x), 1e4, 0, pi)) < 1e-14
    @test abs(radau_quadrature(x -> cos.(x), 1e4, 0, pi)) < 1e-14
    @test abs(chebyshev_quadrature(x -> cos.(x), 1e4, 1, 0, pi)) < 1e-14
    @test abs(chebyshev_quadrature(x -> cos.(x), 1e4, 2, 0, pi)) < 1e-14
end
