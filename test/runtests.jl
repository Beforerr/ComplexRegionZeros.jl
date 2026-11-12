using ComplexRegionZeros
using ComplexRegionZeros: aaa, poles, muller
using Test

@testset "aaa recovers a rational's poles at its exact degree" begin
    ps = [0.21 + 0.1im, -0.5im, -0.7 + 0.05im, 0.9 - 0.15im]
    f(x) = 3 + sum(c / (x - p) for (c, p) in zip((1, 2im, -0.5, 1.0e-3), ps))
    # 51 samples: an odd count; 400: the fit must stop far short of its sample budget
    for z in (complex.(range(-1, 1; length = 51)), [complex(x, y) for x in range(-1, 1; length = 20) for y in range(-0.2, 0.2; length = 20)])
        (; fit, converged) = aaa(f.(z), z)
        @test converged
        @test length(fit.nodes) == length(ps) + 1
        @test fit.(fit.nodes) == f.(fit.nodes)
        @test fit(0.3 + 0.05im) ≈ f(0.3 + 0.05im) rtol = 1.0e-10
        found = poles(fit)
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

@testset "findzeros over a parameter" begin
    f(p) = z -> (z - p) * (z + 1) / (z - 2)
    ps = [0.2, 0.5im, -0.3 + 0.4im]
    rs = findzeros(f, range(-1.5, 1.5, 20) .+ im .* range(-1, 1, 16)', ps, AAA())
    @test all(((r, p),) -> r.converged && minimum(abs, r.roots .- p) < 1.0e-8, zip(rs, ps))
end

@testset "muller" begin
    f(z) = (z - (1 + 2im)) * (z + 3) * exp(z)
    (; root, nevals, converged) = muller(f, 0.8 + 1.7im)
    @test converged && root ≈ 1 + 2im && nevals < 20

    # Non-finite f fails; a trial step into overflow contracts back toward the finite iterates
    @test isnan(muller(x -> ComplexF64(Inf, Inf), 1.0, 1.1, 1.2).root)
    @test isnan(muller(x -> ComplexF64(NaN, NaN), 1.0, 1.1, 1.2).root)
    hits = Ref(0)
    fovf(x) = (v = exp(10x) - 1; isfinite(v) || (hits[] += 1); v)
    @test abs(muller(fovf, -2.0, -1.9, -1.8).root) < 1.0e-10 && hits[] > 0

    # Zero-free, with decades of |f| across any seed triple: a collapsed step is no root.
    slope(x) = exp(1.0e7 * x^2)
    @test !any(muller(slope, 0.002 * cis(θ)).converged for θ in 0:0.2:6.2)

    # A root resolved only to the function's own noise, above `atol`, still converges.
    r0 = 0.3 - 0.1im
    noisy(x) = (x - r0) * (1 + x) + 1.0e-9 * cis(1.0e15 * real(x))
    for seed in (r0, r0 + 0.01)
        r = muller(noisy, seed)
        @test r.converged && abs(r.root - r0) < 1.0e-8
    end
end
