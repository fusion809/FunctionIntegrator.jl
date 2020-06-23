using Integration, Test

# Gaussian curve test
N_gauss_test = 10001;
function expnx2(x)
    exp.(-x.^2)
end

@testset "Gaussian curve" begin
    @test simpsons(expnx2, N_gauss_test, 0, 100) ≈ sqrt(pi)/2
    @test trapezoidal(expnx2, N_gauss_test, 0, 100) ≈ sqrt(pi)/2
    @test legendre_quadrature(expnx2, N_gauss_test, 0, 100) ≈ sqrt(pi)/2
    @test chebyshev_quadrature(expnx2, 1, N_gauss_test, 0, 100) ≈ sqrt(pi)/2
    @test chebyshev_quadrature(expnx2, 2, N_gauss_test, 0, 100) ≈ sqrt(pi)/2
end