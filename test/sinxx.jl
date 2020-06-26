# sin(x)/x function
function sinxx(x)
    return sinc(x/pi);
end

sol_8 = 1.562225466889056293352345138804502677227824980541083456384;

println("Integrating sin(x)/x from 0 to 100 and comparing it to the exact result.")
@testset "sinxx" begin
    println("Running: chebyshev_quadrature with k=1")
    @time @test chebyshev_quadrature(sinxx, 29645, 1, 0, 100) ≈ sol_8
    println("Running: chebyshev_quadrature with k=2")
    @time @test chebyshev_quadrature(sinxx, 41923, 2, 0, 100) ≈ sol_8
    println("Running: chebyshev_quadrature with k=3")
    @time @test chebyshev_quadrature(sinxx, 42083, 3, 0, 100) ≈ sol_8
    println("Running: chebyshev_quadrature with k=4")
    @time @test chebyshev_quadrature(sinxx, 29870, 4, 0, 100) ≈ sol_8
    println("Running: jacobi_quadrature with α=β=1")
    @time @test jacobi_quadrature(sinxx, 65375, 1, 1, 0, 100) ≈ sol_8
    println("Running: legendre_quadrature")
    @time @test legendre_quadrature(sinxx, 37, 0, 100) ≈ sol_8
    println("Running: lobatto_quadrature")
    @time @test lobatto_quadrature(sinxx, 38, 0, 100) ≈ sol_8
    println("Running: radau_quadrature")
    @time @test radau_quadrature(sinxx, 38, 0, 100) ≈ sol_8
    println("Running: rectangle_rule")
    @time @test abs(rectangle_rule(sinxx, 1e8, 0, 100) - sol_8) < 1e-6
    println("Running: simpsons_rule")
    @time @test simpsons_rule(sinxx, 678, 0, 100) ≈ sol_8
    println("Running: trapezoidal_rule")
    @time @test trapezoidal_rule(sinxx, 17622, 0, 100) ≈ sol_8
end
