module ComplexRegionZeros

using LinearAlgebra

export AAA, findzeros
@static if VERSION >= v"1.11.0-DEV.469"
    eval(Meta.parse("public aaa, poles, Barycentric, muller"))
end

include("aaa.jl")
include("muller.jl")

end
