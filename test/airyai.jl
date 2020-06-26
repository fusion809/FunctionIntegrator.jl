# AiryAi(x)
println("Performing the Airy Ai(x) test, where Ai(x) is integrated on the semi-infinite domain, or an approximation of it, namely [0,100] and the result is compared to the analytical result 1/3.")
@testset "Airy Ai(x)" begin
    println("Running: chebyshev_quadrature with k=1")
    @time @test chebyshev_quadrature(x -> airyai(x), 38337, 1, 0, 100) ≈ 1.0/3.0
    println("Running: chebyshev_quadrature with k=2")
    @time @test chebyshev_quadrature(x -> airyai(x), 54215, 2, 0, 100) ≈ 1.0/3.0
    println("Running: chebyshev_quadrature with k=3")
    @time @test chebyshev_quadrature(x -> airyai(x), 54216, 3, 0, 100) ≈ 1.0/3.0
    println("Running: chebyshev_quadrature with k=4")
    @time @test chebyshev_quadrature(x -> airyai(x), 38336, 4, 0, 100) ≈ 1.0/3.0
    println("Running: jacobi_quadrature with α=β=1")
    @time @test jacobi_quadrature(x -> airyai(x), 84543, 1, 1, 0, 100) ≈ 1.0/3.0
    println("Running: laguerre_quadrature with k=1")
    @time @test laguerre_quadrature(x -> airyai(x), 15, 1) ≈ 1.0/3.0
    println("Running: legendre_quadrature")
    @time @test legendre_quadrature(x -> airyai(x), 36, 0, 100) ≈ 1.0/3.0
    println("Running: lobatto_quadrature")
    @time @test lobatto_quadrature(x -> airyai(x), 37, 0, 100) ≈ 1.0/3.0
    println("Running: radau_quadrature")
    @time @test radau_quadrature(x -> airyai(x), 37, 0, 100) ≈ 1.0/3.0
    println("Running: rectangle_rule. Only a rough approximation can be practically achieved using this function.")
    @time @test abs(rectangle_rule(x -> airyai(x), 1e8, 0, 100) - 1.0/3.0) < 1e-6
    println("Running: simpsons_rule")
    @time @test simpsons_rule(x -> airyai(x), 2511, 0, 100) ≈ 1.0/3.0
    println("Running: trapezoidal_rule")
    @time @test trapezoidal_rule(x -> airyai(x), 208381, 0, 100) ≈ 1.0/3.0
end
