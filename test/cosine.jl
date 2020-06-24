# cosine test
@testset "Cosine" begin
    @test abs(rectangle(x -> cos.(x), 1000, 0, pi)) < 1e-14
    @test abs(simpsons(x -> cos.(x), 1000, 0, pi)) < 1e-14
    @test abs(trapezoidal(x -> cos.(x), 1000, 0, pi)) < 1e-14
    @test abs(legendre_quadrature(x -> cos.(x), 1000, 0, pi)) < 1e-14
    @test abs(chebyshev_quadrature(x -> cos.(x), 1, 1000, 0, pi)) < 1e-14
    @test abs(chebyshev_quadrature(x -> cos.(x), 2, 1000, 0, pi)) < 1e-14
end
