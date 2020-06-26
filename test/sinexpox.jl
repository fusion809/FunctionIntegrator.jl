# This is the solution that Wolfram Alpha gave for 0 to inf
sol_11 = 0.3592720446932569240312569249695917230704839918143498570462599122401190492471958577982220255972031016;

function sinexpox(x)
    sin.(x.^2)*exp.(-x).*x^(-1)
end

@testset "sinexpox" begin
    @test chebyshev_quadrature(sinexpox, 655, 1, 0, 100) ≈ sol_11
    @test chebyshev_quadrature(sinexpox, 654, 2, 0, 100) ≈ sol_11
    @test chebyshev_quadrature(sinexpox, 666, 3, 0, 100) ≈ sol_11
    @test chebyshev_quadrature(sinexpox, 576, 4, 0, 100) ≈ sol_11
    @test jacobi_quadrature(sinexpox, 598, 1, 1, 0, 100) ≈ sol_11
    @test laguerre_quadrature(x -> sin.(x.^2).*x.^(-1), 4506, 2) ≈ sol_11
    @test legendre_quadrature(sinexpox, 611, 0, 100) ≈ sol_11
    @test lobatto_quadrature(sinexpox, 598, 1e-14, 100) ≈ sol_11
    @test radau_quadrature(sinexpox, 498, 1e-14, 100) ≈ sol_11
    @test rectangle_rule(sinexpox, 394538, 1e-14, 100) ≈ sol_11
    @test simpsons_rule(sinexpox, 4201, 1e-14, 100) ≈ sol_11
    @test trapezoidal_rule(sinexpox, 394538, 1e-14, 100) ≈ sol_11
end
