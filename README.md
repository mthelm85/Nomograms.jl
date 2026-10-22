# Nomograms

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://mthelm85.github.io/Nomograms.jl/stable)
[![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://mthelm85.github.io/Nomograms.jl/dev)
[![Test workflow status](https://github.com/mthelm85/Nomograms.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/mthelm85/Nomograms.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mthelm85/Nomograms.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mthelm85/Nomograms.jl)
[![Lint workflow Status](https://github.com/mthelm85/Nomograms.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/mthelm85/Nomograms.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/mthelm85/Nomograms.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/mthelm85/Nomograms.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

Build, solve and render **nomograms** — alignment charts that turn an equation into a
graphical calculator, where any straight line laid across the scales picks out a solution.

```julia
using Nomograms

# 1/R = 1/R1 + 1/R2 is already additive: the transform is just 1/x.
nom = additive_nomogram(
    Variable(:R1, ReciprocalScale(); domain = (1, 10), label = "R₁ (Ω)"),
    Variable(:R2, ReciprocalScale(); domain = (1, 10), label = "R₂ (Ω)"),
    Variable(:R,  ReciprocalScale(); label = "R (Ω)");
    title = "Resistors in parallel",
)

evaluate(nom; R1 = 4, R2 = 6)     # (R1 = 4.0, R2 = 6.0, R = 2.4)
evaluate(nom; R1 = 4, R = 2.4)    # ...and backwards: R2 = 6.0

save_svg("parallel.svg", nom; isopleths = [isopleth(nom; R1 = 4, R2 = 6)])
```

Note that `R` declares no domain. The length and position of the middle axis are
*consequences* of the outer two, not free parameters, so the package derives `R ∈ [0.5, 5]`
itself — and rejects a domain you supply that contradicts the relation.

## What it does

- **One geometric core.** Every scale is a parametric curve `value -> (x, y)`, and the
  only thing a nomogram asserts is that one point per scale is collinear (d'Ocagne's
  determinant). Solving is therefore always "where does this line cross that curve",
  which means forward and inverse problems are the same call and new layouts need no new
  solver.
- **No silent wrong answers on curved scales.** A straight line can cross a curve more
  than once, and each crossing is a real solution. The whole domain is searched rather
  than the residual's sign being checked at the two ends, which would report *no* solution
  for an even number of crossings. `evaluate` refuses to pick one of several and names
  them; `evaluate_all` returns the lot. Straight scales have a monotone residual, so
  parallel layouts never see this.
- **Invertible scale transforms.** `LinearScale`, `LogScale`, `ReciprocalScale`,
  `PowerScale`, `ComposedScale`, and `FunctionScale(f, finv)` for anything else. Carrying
  the inverse alongside the forward map is what lets the result axis be derived exactly.
- **Readable axes.** Ticks are chosen as round numbers *in the variable* and then thinned
  by distance *on the page*, so a reciprocal or logarithmic axis keeps sensible labels
  where it is stretched and quietly drops them where it is compressed. Endpoints and
  decade boundaries are kept in preference to crowded neighbours.
- **Two renderers over one layout.** `drawing(nom)` produces pure geometry; `svg` walks it
  with no dependencies beyond the standard library, and a package extension walks it with
  Makie when you load a backend.

## Rendering

The built-in SVG writer needs nothing extra:

```julia
save_svg("chart.svg", nom; isopleths = [...], pxperunit = 4, background = "none")
```

Load any Makie backend for a `Figure`, and with it PDF, PNG and Makie's own `save`:

```julia
using CairoMakie
fig = nomogramplot(nom; isopleths = [isopleth(nom; R1 = 4, R2 = 6)])
save("chart.pdf", fig)
```

Both go through `drawing(nom; style = NomogramStyle(...))`, so a `NomogramStyle` means the
same thing in either backend, and tests can assert on tick positions without rendering.

## Products, and anything separable

A product relation becomes additive under a logarithm, which
`multiplicative_nomogram` applies for you:

```julia
nom = multiplicative_nomogram(
    Variable(:rate;  domain = (1, 100)),
    Variable(:hours; domain = (1, 100)),
    Variable(:cost),                        # derived: 1 to 10 000
)
evaluate(nom; rate = 7.5, hours = 12).cost  # 90.0
```

More generally, a three-variable relation is chartable on parallel scales exactly when it
can be written `f₁(u) + f₂(v) = f₃(w)`. Supply the `fᵢ` as transforms and the construction
does the rest; see [`examples/compound_growth.jl`](examples/compound_growth.jl), which
charts `growth = (1 + r)ⁿ` by way of `log(growth) = n · log(1 + r)`.

For a layout the constructions do not cover, build `Scale`s directly — a curve per
variable, straight or not — and pass them to `Nomogram`. Everything downstream, the solver
included, works unchanged.

## Examples

- [`examples/parallel_resistance.jl`](examples/parallel_resistance.jl) — resistors in
  parallel, the classic reciprocal chart.
- [`examples/compound_growth.jl`](examples/compound_growth.jl) — compound growth, with a
  custom transform and hand-picked ticks for an axis that compresses hard.

## Status

Version 0.1.0, not yet registered. The parallel-axis construction, the solver, both
renderers and the tick engine are implemented and tested. Not yet implemented: non-parallel
layouts (N/Z charts, curved and circular scales) as named constructions, grid or
"compound" nomograms for four or more variables, and a symbolic front end that would
recognise separable forms from an expression rather than taking the transforms directly.
