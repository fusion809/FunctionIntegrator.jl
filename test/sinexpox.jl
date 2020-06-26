# This is the solution that Wolfram Alpha gave for 0 to inf
sol_11 = 0.3592720446932569240312569249695917230704839918143498570462599122401190492471958577982220255972031016;

function sinexpox(x)
    if x!=0
        sin.(x.^2)*exp.(-x).*x^(-1)
    else
        0
    end
end

println("Integrating sin(x^2)e^(-x)/x from 0 to infinity, with the approximated domain of integration of 0 to 100. Integrand function has been adjusted so that the singularity at x=0 no longer exists.")
@testset "sinexpox" begin
    println("Running: chebyshev_quadrature with k=1")
    @time @test chebyshev_quadrature(sinexpox, 655, 1, 0, 100) ≈ sol_11
    println("Running: chebyshev_quadrature with k=2")
    @time @test chebyshev_quadrature(sinexpox, 654, 2, 0, 100) ≈ sol_11
    println("Running: chebyshev_quadrature with k=3")
    @time @test chebyshev_quadrature(sinexpox, 666, 3, 0, 100) ≈ sol_11
    println("Running: chebyshev_quadrature with k=4")
    @time @test chebyshev_quadrature(sinexpox, 576, 4, 0, 100) ≈ sol_11
    println("Running: jacobi_quadrature with α=β=1")
    @time @test jacobi_quadrature(sinexpox, 598, 1, 1, 0, 100) ≈ sol_11
    println("Running: laguerre_quadrature with k=2")
    @time @test laguerre_quadrature(x -> sin.(x.^2).*x.^(-1), 4506, 2) ≈ sol_11
    println("Running: legendre_quadrature")
    @time @test legendre_quadrature(sinexpox, 611, 0, 100) ≈ sol_11
    println("Running: lobatto_quadrature")
    @time @test lobatto_quadrature(sinexpox, 598, 0, 100) ≈ sol_11
    println("Running: radau_quadrature")
    @time @test radau_quadrature(sinexpox, 498, 0, 100) ≈ sol_11
    println("Running: rectangle_rule")
    @time @test rectangle_rule(sinexpox, 394538, 0, 100) ≈ sol_11
    println("Running: simpsons_rule")
    @time @test simpsons_rule(sinexpox, 4201, 0, 100) ≈ sol_11
    println("Running: trapezoidal_rule")
    @time @test trapezoidal_rule(sinexpox, 394538, 0, 100) ≈ sol_11
end
