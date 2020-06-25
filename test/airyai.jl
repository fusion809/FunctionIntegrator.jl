# AiryAi(x)
@testset "Airy Ai(x)" begin
    @test abs(rectangle_rule(x -> airyai(x), 1e8, 0, 100) - 1.0/3.0) < 1e-6
    @test simpsons_rule(x -> airyai(x), 2511, 0, 100) ≈ 1.0/3.0
    @test trapezoidal_rule(x -> airyai(x), 208381, 0, 100) ≈ 1.0/3.0
    @test jacobi_quadrature(x -> airyai(x), 84543, 1, 1, 0, 100) ≈ 1.0/3.0
    @test laguerre_quadrature(x -> airyai(x), 15, 1) ≈ 1.0/3.0
    @test legendre_quadrature(x -> airyai(x), 36, 0, 100) ≈ 1.0/3.0
    @test lobatto_quadrature(x -> airyai(x), 37, 0, 100) ≈ 1.0/3.0
    @test radau_quadrature(x -> airyai(x), 37, 0, 100) ≈ 1.0/3.0
    @test chebyshev_quadrature(x -> airyai(x), 38337, 1, 0, 100) ≈ 1.0/3.0
    @test chebyshev_quadrature(x -> airyai(x), 54215, 2, 0, 100) ≈ 1.0/3.0
end
