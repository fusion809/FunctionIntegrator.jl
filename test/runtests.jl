using Integration
using Test

function f(x)
    exp.(-x.^2)
end

@testset "general" begin
    @test simpsons(f, 1001, 0, 100) ≈ sqrt(pi)/2
end