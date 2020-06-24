using SpecialFunctions;

sol_9 = exp(1)*pi*erfc(1);
@testset "expnx2datan" begin
    @test rectangle(x -> exp.(-x.^2).*(x.^2+1).^(-1), 1e8, -100, 100) ≈ sol_9
    @test simpsons(x -> exp.(-x.^2).*(x.^2+1).^(-1), 10000, -100, 100) ≈ sol_9
    @test trapezoidal(x -> exp.(-x.^2).*(x.^2+1).^(-1), 10000, -100, 100) ≈ sol_9
    @test hermite_quadrature(x -> (x.^2+1).^(-1), 100, 2) ≈ sol_9
    @test laguerre_quadrature(x -> 2*exp.(-x.^2).*(x.^2+1).^(-1), 100, 1) ≈ sol_9
    @test legendre_quadrature(x -> exp.(-x.^2).*(x.^2+1).^(-1), 10000, -100, 100) ≈ sol_9
    @test chebyshev_quadrature(x -> exp.(-x.^2).*(x.^2+1).^(-1), 1, 100000, -100, 100) ≈ sol_9
    @test chebyshev_quadrature(x -> exp.(-x.^2).*(x.^2+1).^(-1), 2, 100000, -100, 100) ≈ sol_9
end
