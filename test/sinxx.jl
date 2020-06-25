# sin(x)/x function
function sinxx(x)
    return sinc(x/pi);
end

sol_8 = 1.562225466889056293352345138804502677227824980541083456384;

@testset "sinxx" begin
    @test abs(rectangle(sinxx, 1e8, 0, 100) - sol_8) < 1e-6
    @test simpsons(sinxx, 1e5, 0, 100) ≈ sol_8
    @test trapezoidal(sinxx, 1e5, 0, 100) ≈ sol_8
    @test jacobi_quadrature(sinxx, 1e6, 1, 1, 0, 100) ≈ sol_8
    @test legendre_quadrature(sinxx, 1e6, 0, 100) ≈ sol_8
    @test lobatto_quadrature(sinxx, 1e6, 0, 100) ≈ sol_8
    @test radau_quadrature(sinxx, 1e6, 0, 100) ≈ sol_8
    @test chebyshev_quadrature(sinxx, 1e5, 1, 0, 100) ≈ sol_8
    @test chebyshev_quadrature(sinxx, 1e5, 2, 0, 100) ≈ sol_8
end
