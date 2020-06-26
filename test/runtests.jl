# These tests are run in Travis CI:
# https://travis-ci.com/github/fusion809/FunctionIntegrator.jl
using FunctionIntegrator, SpecialFunctions, Test

include("gaussian.jl")
include("simppen.jl")
include("cosine.jl")
include("logarithm.jl")
include("cos_cot_integral.jl")
include("logoverx.jl")
include("test_7.jl")
include("sinxx.jl")
include("airyai.jl")
include("expnx2datan.jl")
include("modbessel0.jl")
include("sinexpox.jl")
include("besselj.jl")
