# sin(x)/x function
function sinxx(x)
    return sinc(x/pi);
end

sol_8 = 1.562225466889056293352345138804502677227824980541083456384;

@testset "sinxx" begin
    @test chebyshev_quadrature(sinxx, 29645, 1, 0, 100) ≈ sol_8
    @test chebyshev_quadrature(sinxx, 41923, 2, 0, 100) ≈ sol_8
    @test chebyshev_quadrature(sinxx, 42083, 3, 0, 100) ≈ sol_8
    @test chebyshev_quadrature(sinxx, 29870, 4, 0, 100) ≈ sol_8
    @test jacobi_quadrature(sinxx, 65375, 1, 1, 0, 100) ≈ sol_8
    @test legendre_quadrature(sinxx, 37, 0, 100) ≈ sol_8
    @test lobatto_quadrature(sinxx, 38, 0, 100) ≈ sol_8
    @test radau_quadrature(sinxx, 38, 0, 100) ≈ sol_8
    @test abs(rectangle_rule(sinxx, 1e8, 0, 100) - sol_8) < 1e-6
    @test simpsons_rule(sinxx, 678, 0, 100) ≈ sol_8
    @test trapezoidal_rule(sinxx, 17622, 0, 100) ≈ sol_8
end
