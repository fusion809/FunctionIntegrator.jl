# AiryAi(x)
@testset "Airy Ai(x)" begin
    @test abs(rectangle_rule(x -> airyai(x), 1e8, 0, 100) - 1.0/3.0) < 1e-6
    @test simpsons_rule(x -> airyai(x), 1e6, 0, 100) ≈ 1.0/3.0
    @test trapezoidal_rule(x -> airyai(x), 1e6, 0, 100) ≈ 1.0/3.0
    @test abs(jacobi_quadrature(x -> airyai(x), 1e3, 1, 1, 0, 100) - 1.0/3.0) < 1e-3
    @test laguerre_quadrature(x -> airyai(x), 120, 1) ≈ 1.0/3.0
    @test legendre_quadrature(x -> airyai(x), 1e3, 0, 100) ≈ 1.0/3.0
    @test lobatto_quadrature(x -> airyai(x), 1e3, 0, 100) ≈ 1.0/3.0
    @test radau_quadrature(x -> airyai(x), 1e3, 0, 100) ≈ 1.0/3.0
    @test chebyshev_quadrature(x -> airyai(x), 1e5, 1, 0, 100) ≈ 1.0/3.0
    @test chebyshev_quadrature(x -> airyai(x), 1e5, 2, 0, 100) ≈ 1.0/3.0
end
