function muller(f, x0; kw...)
    h = _seed_offset(x0)
    return muller(f, x0 - h, x0 + h, x0 + h * im; kw...)
end

# Absolute-floored relative step capped at 0.1|x0|; x0 = 0 has no scale ⇒ a small absolute step.
_seed_offset(x0) = iszero(x0) ? 1.0e-3 : min(1.0e-3 * max(abs(x0), 1.0), 0.1 * abs(x0))

function muller(f, x0, x1, x2; atol = 1.0e-10, maxiter = 100)
    x0, x1, x2 = promote(complex(x0), complex(x1), complex(x2))
    nan = typeof(x2)(NaN, NaN)
    f0, f1, f2 = f(x0), f(x1), f(x2)
    a0, a1, a2 = abs(f0), abs(f1), abs(f2)
    nevals = 3
    fail(nevals) = (; root = nan, nevals, converged = false)
    fell = false
    rtol = sqrt(atol)
    for _ in 1:maxiter
        h1 = x1 - x0
        h2 = x2 - x1
        (h1 == 0 || h2 == 0 || h1 + h2 == 0) && return fail(nevals)
        delta1 = (f1 - f0) / h1
        delta2 = (f2 - f1) / h2
        a = (delta2 - delta1) / (h2 + h1)
        b = a * h2 + delta2
        c = f2

        disc = sqrt(b^2 - 4a * c)
        # The larger of |b ± disc|: |b + disc|² − |b − disc|² = 4 Re(b̄·disc)
        denom = real(conj(b) * disc) > 0 ? b + disc : b - disc
        denom == 0 && return fail(nevals)

        dx = -2c / denom
        x3 = x2 + dx
        f3 = f(x3)
        a3 = abs(f3)
        nevals += 1
        # Backtrack steps whose value is non-finite or blows up past the iterates'.
        flimit = 1.0e6 * max(a0, a1, a2)
        for _ in 1:12
            (isfinite(f3) && a3 <= flimit) && break
            dx /= 2
            x3 = x2 + dx
            f3 = f(x3)
            a3 = abs(f3)
            nevals += 1
        end
        isfinite(f3) || return fail(nevals)
        # A small step alone is no root: exponential decay mimics one, widely differing |f| collapse
        # the step without improving f, and absolute residuals depend on f's normalization. So also
        # require |f| to have fallen substantially — and remember that it did, since evaluation
        # noise near a root can stop further decrease.
        mid = max(min(a0, a1), min(max(a0, a1), a2))
        step = abs(dx) / max(abs(x3), 1)
        fell = (a3 <= a2 / 10 && a3 <= rtol * mid) || (fell && step <= rtol)
        (iszero(f3) || (fell && step <= atol)) && return (; root = x3, nevals, converged = true)

        x0, x1, x2 = x1, x2, x3
        f0, f1, f2 = f1, f2, f3
        a0, a1, a2 = a1, a2, a3
    end
    return fail(nevals)
end
