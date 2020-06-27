# These tests are run in Travis CI:
# https://travis-ci.com/github/fusion809/FunctionIntegrator.jl
using FunctionIntegrator, SpecialFunctions, Test

include("airyai.jl")
include("besselj.jl")
include("cos_cot_integral.jl")
include("cosine.jl")
include("expnx2datan.jl")
include("gaussian.jl")
include("logarithm.jl")
include("logoverx.jl")
include("modbessel0.jl")
include("simppen.jl")
include("sinexpox.jl")
include("sinxx.jl")
include("test_7.jl")
