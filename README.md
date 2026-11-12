# ComplexRegionZeros.jl

All zeros of a meromorphic function inside a region of the complex plane, without starting guesses.

```julia
using ComplexRegionZeros

f(z) = (z^2 - 1) * (z - 0.5im) / (z - 2)
(; roots, nevals, converged) = findzeros(f, range(-3, 3, 20) .+ im .* range(-1, 1, 16)', AAA())
```

`roots` are candidates, whatever the algorithm: they may include spurious ones (Froissart
doublets) and ones outside the sampled region, left for the caller to filter by its own region
and residual. Each can be refined by a local method such as Muller or Newton.

## API

Exported: `findzeros`, `AAA`. Public but not exported: `aaa`, `Barycentric`, `poles`, `muller`.

- `findzeros(f, domain, alg) -> (; roots, nevals, converged)`: zeros of `f` over `domain`, in
  the form `alg` consumes; `AAA` takes an array of sample points. `converged` is `false` when
  the fit stopped short of `tol` and may be missing zeros, including when no sample is usable.
- `AAA(; tol = 1e-13, max_degree = 150, stagnation = 5)`: fits `1/f` by AAA and returns the
  fit's poles. Non-finite and zero values are dropped, so `f` can return `NaN` where its
  value cannot be trusted; fitting such values would spend the fit's degree on their error.
- `aaa(y, z; tol, max_degree, stagnation) -> (; fit, converged)`: discrete AAA fit of values
  `y` at points `z`. It stops on the same rule as RationalFunctionApproximation's
  `approximate(y, z, Barycentric())`, falling back on the best fit when the error stagnates
  or `max_degree` is hit.
- `Barycentric(nodes, values, weights)`: barycentric rational, callable as `r(z)`.
- `poles(r::Barycentric)`: its finite poles.
- `muller(f, x0; atol = 1e-10, maxiter = 100)`, `muller(f, x0, x1, x2; …) -> (; root, nevals, converged)`:
  Muller's method from a guess or three starting points; `root` is `NaN` unless converged. A
  step below `atol·max(|x|, 1)` converges only once `|f|` has fallen by `√atol`, so exponential
  decay or a collapsed step is not taken for a root; steps to non-finite or overflowing `f` are
  halved back.

### Extending

A new algorithm is a type plus a `findzeros(f, domain, ::Alg)` method; its `domain` may be
a geometric region instead of points, and its result may carry extra fields. Problems over a
parameter (e.g. a wavevector sweep) wrap `findzeros(z -> f(z, p), domain, alg)` per `p` and
link the roots into branches downstream, as VlasovMaxwellDispersion does; a fit over `(z, p)`
jointly would be its own `findzeros` method that reuses `poles` on each `p`-slice.

The fit and its cost are derived in [docs/aaa.md](docs/aaa.md). `benchmark/compare.jl`
compares fits against RationalFunctionApproximation.jl and BaryRational.jl.
