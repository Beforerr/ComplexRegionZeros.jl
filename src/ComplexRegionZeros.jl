module ComplexRegionZeros

using LinearAlgebra

export AAA, findzeros
@static if VERSION >= v"1.11.0-DEV.469"
    eval(Meta.parse("public aaa, poles, Barycentric, muller"))
end

# Algorithms that share work across the parameter dimension specialize this.
findzeros(f, domain, ps, alg) = map(p -> findzeros(f(p), domain, alg), ps)

include("aaa.jl")
include("muller.jl")

end
