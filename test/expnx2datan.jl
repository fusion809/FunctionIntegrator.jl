# Approximating the integral from -inf to inf of e^(-x^2)/(x^2+1)
sol_9 = exp(1)*pi*erfc(1);
@testset "expnx2datan" begin
    @test chebyshev_quadrature(x -> exp.(-x.^2).*(x.^2+1).^(-1), 1029, 1, -100, 100) ≈ sol_9
    @test chebyshev_quadrature(x -> exp.(-x.^2).*(x.^2+1).^(-1), 1028, 2, -100, 100) ≈ sol_9
    @test hermite_quadrature(x -> (x.^2+1).^(-1), 53, 2) ≈ sol_9
    @test jacobi_quadrature(x -> exp.(-x.^2).*(x.^2+1).^(-1), 1027, 1, 1, -100, 100) ≈ sol_9
    @test laguerre_quadrature(x -> 2*exp.(-x.^2).*(x.^2+1).^(-1), 55, 1) ≈ sol_9
    @test legendre_quadrature(x -> exp.(-x.^2).*(x.^2+1).^(-1), 1028, -100, 100) ≈ sol_9
    @test lobatto_quadrature(x -> exp.(-x.^2).*(x.^2+1).^(-1), 1029, -100, 100) ≈ sol_9
    @test radau_quadrature(x -> exp.(-x.^2).*(x.^2+1).^(-1), 669, -100, 100) ≈ sol_9
    @test rectangle_rule(x -> exp.(-x.^2).*(x.^2+1).^(-1), 655, -100, 100) ≈ sol_9
    @test simpsons_rule(x -> exp.(-x.^2).*(x.^2+1).^(-1), 1233, -100, 100) ≈ sol_9
    @test trapezoidal_rule(x -> exp.(-x.^2).*(x.^2+1).^(-1), 655, -100, 100) ≈ sol_9
end
