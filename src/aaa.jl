Base.@kwdef struct AAA
    tol::Float64 = 1.0e-13
    max_degree::Int = 150
    stagnation::Int = 5
end

function findzeros(f, z::AbstractArray, alg::AAA)
    Z = vec(z)
    nevals = length(Z)
    F = map(f, Z)
    nvalid = 0
    @inbounds for v in F
        nvalid += isfinite(v) && !iszero(v)
    end
    if nvalid == length(F)
        map!(inv, F, F)
    else
        z = similar(Z, nvalid)
        y = similar(F, nvalid)
        i = 0
        @inbounds for j in eachindex(F, Z)
            v = F[j]
            (isfinite(v) && !iszero(v)) || continue
            i += 1
            z[i], y[i] = Z[j], inv(v)
        end
        Z, F = z, y
    end
    (; fit, converged) = aaa(F, Z; alg.tol, alg.max_degree, alg.stagnation)
    return (; roots = poles(fit), nevals, converged)
end

struct Barycentric{Z <: AbstractVector, V <: AbstractVector, W <: AbstractVector}
    nodes::Z
    values::V
    weights::W
end

function (r::Barycentric)(x)
    (; nodes, values, weights) = r
    T = promote_type(typeof(x), eltype(nodes), eltype(values), eltype(weights))
    num = den = zero(T)
    for k in eachindex(nodes, values, weights)
        δ = x - nodes[k]
        iszero(δ) && return T(values[k])
        c = weights[k] / δ
        num += c * values[k]
        den += c
    end
    return num / den
end

# Poles of a barycentric rational = finite eigenvalues of the arrowhead pencil (Klein 2012);
# the infinite ones are the pencil's own, not the fit's. Zero weights leave the pencil singular,
# so drop those nodes first.
function poles(r::Barycentric)
    keep = map(!iszero, r.weights)
    z, w = r.nodes[keep], r.weights[keep]
    m = length(w)
    T = real(eltype(z))
    B = diagm([zero(T); ones(T, m)])
    E = [zero(T) transpose(w); ones(T, m) diagm(z)]
    return filter!(isfinite, eigvals(E, B))
end

# The Loewner matrix changes by one column and one zeroed row (the new node's) per step, so
# its thin QR is updated rather than refactored and the weights come from inverse iteration
# on R: O(m·n) per step where a fresh SVD costs O(m·n²). Derivation: docs/aaa.md.
function aaa(y, z; tol = 1.0e-13, max_degree = 150, stagnation = 5)
    m = length(z)
    T = promote_type(eltype(z), eltype(y))
    m == 0 && return (; fit = Barycentric(z[1:0], y[1:0], T[]), converged = false)
    fmax = maximum(abs, y)
    # A node's own Loewner row is dropped, so the fit can never outrank the surviving samples.
    N = min(max_degree + 1, m >> 1)
    C = Matrix{T}(undef, m, N)      # Cauchy 1/(zᵢ − σₖ)
    Q = Matrix{T}(undef, m, N)      # Loewner (yᵢ − fσₖ)·C over the test rows = Q·R
    R = zeros(T, N, N)
    a = Vector{T}(undef, m); u = Vector{T}(undef, m); h = Vector{T}(undef, N)
    num = Vector{T}(undef, m); den = Vector{T}(undef, m)
    test = trues(m)
    _, i0 = findmin(abs, y)         # first node: same seed as `approximate`
    σ = [i0]; test[i0] = false
    errs = Float64[]
    ws = Vector{Vector{T}}()        # weights of every step; step s interpolates at σ[1:s]
    for n in 1:N
        j = σ[n]
        _qr_zero_row!(Q, R, j, n - 1, u, h)
        @inbounds for i in 1:m
            δ = z[i] - z[j]
            C[i, n] = iszero(δ) ? 1 / eps() : 1 / δ   # δ=0 only on duplicate samples
            a[i] = test[i] ? (y[i] - y[j]) * C[i, n] : zero(T)
        end
        _qr_append!(Q, R, a, n - 1, h)
        w = _nullvector!(Vector{T}(undef, n), R)

        Cn = view(C, :, 1:n)
        hn = view(h, 1:n)
        hn .= w .* view(y, σ)
        mul!(num, Cn, hn)
        mul!(den, Cn, w)
        err, k = -Inf, 0
        @inbounds for i in 1:m
            test[i] || continue
            e = abs(y[i] - num[i] / den[i])
            e > err && ((err, k) = (e, i))
        end
        push!(errs, err)
        push!(ws, w)

        s = _quitting_check(errs, stagnation, tol, fmax, N)
        # s > 0, unconverged: fall back on the best fit
        if s != 0
            b = s < 0 ? n : s
            return (; fit = Barycentric(z[σ[1:b]], y[σ[1:b]], ws[b]), converged = s < 0)
        end
        push!(σ, k); test[k] = false
    end
    # A single sample: nothing to fit.
    return (; fit = Barycentric(z[σ], y[σ], ones(T, length(σ))), converged = false)
end

# Thin QR `A = Q[:, 1:n]·R[1:n, 1:n]`, updated in place as `A` gains a column or has a row
# zeroed. `h` (length ≥ n) and `u` (length m) are scratch.

# Append column `a` (destroyed): classical Gram–Schmidt, repeated while the remainder
# keeps collapsing — one pass leaves it unorthogonal exactly when `a` is nearly dependent.
function _qr_append!(Q, R, a, n, h)
    Qn = view(Q, :, 1:n); r = view(R, 1:n, n + 1); hn = view(h, 1:n)
    fill!(r, 0)
    nrm = norm(a)
    for _ in 1:3
        mul!(hn, Qn', a)
        mul!(a, Qn, hn, -1, 1)
        r .+= hn
        prev, nrm = nrm, norm(a)
        nrm > prev / 2 && break
    end
    R[n + 1, n + 1] = nrm
    q = view(Q, :, n + 1)
    iszero(nrm) ? fill!(q, 0) : (q .= a ./ nrm)
    return
end

# Zero row `p` of `A`: extend Q by the part `u` of e_p outside its range, then rotate columns
# n…1 against `u` until row p of Q vanishes. The same rotations keep R triangular and gather
# the removed row in a discarded extra row of R (`h`).
function _qr_zero_row!(Q, R, p, n, u, h)
    n == 0 && return
    Qn = view(Q, :, 1:n); hn = view(h, 1:n)
    fill!(u, 0); u[p] = 1
    for _ in 1:2
        mul!(hn, Qn', u)
        mul!(u, Qn, hn, -1, 1)
    end
    τ = norm(u)
    iszero(τ) || (u ./= τ)
    fill!(hn, 0)
    b = real(u[p])
    @inbounds for j in n:-1:1
        q = Q[p, j]
        ρ = sqrt(abs2(q) + b^2)
        iszero(ρ) && continue
        c, s = b / ρ, q / ρ
        for i in axes(Q, 1)
            qi, ui = Q[i, j], u[i]
            Q[i, j] = c * qi - s * ui
            u[i] = conj(s) * qi + c * ui
        end
        Q[p, j] = 0
        for l in j:n
            rj, rh = R[j, l], h[l]
            R[j, l] = c * rj - conj(s) * rh
            h[l] = s * rj + c * rh
        end
        b = ρ
    end
    return
end

# Right singular vector of the smallest singular value of the upper-triangular `R`, by
# inverse iteration on RᴴR. A few sweeps reach a vector whose `‖R·w‖` is the smallest
# singular value to within a small factor, which is all the fit asks of its weights; the
# floored diagonal keeps an exactly singular `R` (a reproduced rational) solvable.
function _nullvector!(w, R)
    n = length(w)
    T = real(eltype(R))
    tiny = max(eps(T) * maximum(i -> abs(R[i, i]), 1:n), sqrt(floatmin(T)))
    d(i) = (x = R[i, i]; abs(x) < tiny ? oftype(x, tiny) : x)
    fill!(w, 1)
    @inbounds for _ in 1:4
        for i in 1:n                 # Rᴴx = w
            s = w[i]
            for l in 1:(i - 1)
                s -= conj(R[l, i]) * w[l]
            end
            w[i] = s / conj(d(i))
        end
        normalize!(w)
        for i in n:-1:1              # Rw = x
            s = w[i]
            for l in (i + 1):n
                s -= R[i, l] * w[l]
            end
            w[i] = s / d(i)
        end
        normalize!(w)
    end
    return w
end

# `RationalFunctionApproximation.quitting_check`, reproduced so the fit stops exactly where
# `approximate` did: -1 converged (keep the last fit), 0 continue, n rewind to fit n.
function _quitting_check(errs, stagnation, tol, fmax, max_iter)
    n = length(errs)
    errs[end] <= tol * fmax && return -1
    stagnant = false
    min_err, min_k = findmin(errs)
    if n >= stagnation + 5 && min_err < fmax / 100
        plateau = 5min_err
        stagnant = all(plateau < e for e in last(errs, stagnation)) || (min_k < n - 2stagnation)
    end
    return (n >= max_iter || stagnant) ? argmin(errs) : 0
end
