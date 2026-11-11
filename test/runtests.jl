using ComplexRegionZeros
using ComplexRegionZeros: aaa, poles
using Test

@testset "aaa recovers a rational's poles at its exact degree" begin
    ps = [0.21 + 0.1im, -0.5im, -0.7 + 0.05im, 0.9 - 0.15im]
    f(x) = 3 + sum(c / (x - p) for (c, p) in zip((1, 2im, -0.5, 1.0e-3), ps))
    # 51 samples: an odd count; 400: the fit must stop far short of its sample budget
    for z in (complex.(range(-1, 1; length = 51)), [complex(x, y) for x in range(-1, 1; length = 20) for y in range(-0.2, 0.2; length = 20)])
        r = aaa(f.(z), z)
        @test r.converged
        @test length(r.nodes) == length(ps) + 1
        @test r.values == f.(r.nodes)
        found = poles(r.nodes, r.weights)
        @test all(p -> minimum(abs, found .- p) < 1.0e-9, ps)
    end
end

@testset "findzeros skips non-finite samples" begin
    zs = [1.0, -1.0, 0.5im]
    f(z) = prod(z .- zs) / (z - 2)
    # the 13×5 grid holds z = 2 exactly, where f is infinite
    Z = range(-3, 3, 13) .+ im .* range(-1, 1, 5)'
    r = findzeros(f, Z, AAA())
    @test r.converged
    @test r.nevals == length(Z)
    @test all(z -> minimum(abs, r.roots .- z) < 1.0e-8, zs)
    r = findzeros(z -> Inf, Z, AAA())
    @test isempty(r.roots) && !r.converged
end

@testset "findzeros on a filtered grid and on boundary points" begin
    zs = [1.0, -1.0, 0.5im]
    f(z) = prod(z .- zs) / (z - 2)
    for Z in (filter(z -> abs(z) < 1.5, range(-1.5, 1.5, 30) .+ im .* range(-1.5, 1.5, 30)'), 1.5 .* cispi.(2 .* (0:99) ./ 100))
        r = findzeros(f, Z, AAA())
        @test all(z -> minimum(abs, r.roots .- z) < 1.0e-8, zs)
    end
end
