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

## Or just give it the formula

Every constructor below asks you to rewrite your formula into a separable form first and
choose the transforms. `auto_nomogram` does that for you, from the function as you would
ordinarily write it:

```julia
nom = auto_nomogram((r, n) -> (1 + r/100)^n, :r => (1, 20), :n => (1, 30), :growth)

evaluate(nom; r = 5, n = 10).growth   # 1.6288946…  = 1.05^10
evaluate(nom; r = 7, growth = 2).n    # 10.2447…    the rule of 72's question
```

It applies the classical test for when a formula `z = F(x, y)` can be drawn on three
straight scales at all — transforms with `h(z) = f(x) + g(y)` exist exactly when
`ln|F_x / F_y|` splits into a part in `x` plus a part in `y` — and when it can, the same
identity hands back the transforms by integration. Each axis then gets linear or
decade-based ticks according to how its transform spreads values out.

Nothing is returned on trust. The finished chart is *solved* over a grid of inputs and
every answer compared with your function; a chart that disagrees is an error, not a
result. When a formula has no straight-scale chart, it says so and points you at the
constructions that can still handle it.

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
- **Labels that stay out of each other's way.** Each label's text box is computed from its
  anchor, alignment and rotation, so the layout can see collisions between *different*
  scales — the ones that per-axis thinning cannot know about. Axis titles win, tick labels
  that would overlap are dropped while their tick marks stay, and isopleth labels are
  nudged clear rather than discarded. Titles follow their scale off its end, and
  `ticklabelrotation = :tangent` turns labels to run along a steep or curved axis.
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

The defaults follow the proportions of a classic d3 axis — hairline rules, short ticks,
small type. All lengths are in chart units, and at the default `pxperunit = 4` they come
out as 1px strokes, 6px major ticks and 10px labels.

Isopleths get **three** colours, not more: blue, orange, aqua. Every isopleth on a chart is
visible at once, so all three pairs must stay apart under colour-vision deficiency, and no
fourth hue holds up once every pair is in play. These were checked with a validator rather
than by eye — worst pair ΔE 9.2 under deuteranopia, 24.0 in normal vision, against a ≥8
target. d3's own `category10` is *not* used for a measurable reason: its orange and green
sit ΔE 0.7 apart under protanopia, which is to say they are the same colour. Past three
isopleths the dash pattern advances instead of the hue, so identity never rests on colour
alone.

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

## Ratios, on an N-chart

Parallel scales are not the only layout. `ratio_nomogram` builds an N-chart (Z-chart) —
two vertical scales joined by a diagonal, with the result read off the diagonal — for

```julia
nom = ratio_nomogram(
    Variable(:distance; domain = (1, 10)),
    Variable(:time;     domain = (1, 10)),
    Variable(:speed),                        # derived: 0.1 to 10
)
evaluate(nom; distance = 6, time = 3).speed  # 2.0
```

`multiplicative_nomogram` charts ratios too, by taking logarithms. The N-chart earns its
place by needing none: the outer scales stay **linear**, so precision is even along them
and small ratios are not squeezed into a corner.

The diagonal is a projective scale, `s = k·r / (k·r + 1)`, so a ratio spanning decades is
spread across it roughly logarithmically — and the tick values follow suit rather than
being picked linearly, which over 0.1 to 10 would offer nothing between 0.1 and 2.

## Anything else, from the determinant

`determinant_nomogram` takes the scales as curves and fits them to the canvas, which covers
the layouts the named constructions cannot express. A quadratic needs a curved scale:

```julia
# Matching (1-x₃)·p + x₃·q = y₃ against z·p + q + z² = 0 gives the curve below.
nom = determinant_nomogram(
    Variable(:p; domain = (-6, 0)) => p -> (0.0, p),
    Variable(:q; domain = (0, 6))  => q -> (1.0, q),
    Variable(:z; domain = (0.2, 5)) => z -> (1 / (1 + z), -z^2 / (1 + z));
    relation = (p, q, z) -> z^2 + p * z + q,
)

[sol.z for sol in evaluate_all(nom; p = -3, q = 2)]   # [1.0, 2.0]
```

Write the curves in whatever coordinates suit the geometry; their bounding box is measured
and mapped onto the canvas. Collinearity survives any affine map, so the fit cannot change
what the chart computes.

Note the two answers. A quadratic has two roots and the chart shows why — one alignment
line, meeting the curved scale twice.

Getting a determinant form right is the hard part of designing a chart by hand, and a wrong
one still draws perfectly well, so `max_relation_residual(nom)` solves on the geometry and
checks the answers against the algebra.

For full control, build `Scale`s directly — a curve per variable, in final coordinates —
and pass them to `Nomogram`. Everything downstream, the solver included, works unchanged.

## More than three variables

One alignment line can relate only three scales, so a four-variable relation is split into
two stages sharing an intermediate, and `compound` positions them so the shared scale
coincides:

```julia
# V = L × W × H, as t = L × W then V = t × H
stage1 = multiplicative_nomogram(L, W, Variable(:t; ticks = Float64[]))
stage2 = multiplicative_nomogram(
    Variable(:t; domain = stage1[:t].domain, ticks = Float64[]), H, V,
)
nom = compound(stage1, stage2)

evaluate(nom; L = 3, W = 4, H = 2)     # (L=3, t=12, W=4, V=24, H=2)
evaluate(nom; L = 3, W = 4, V = 24).H  # 2.0 — any direction the chain allows
```

The reader lays a line across the first stage, notes where it crosses the shared scale,
then pivots about that point and reads the second. Hence the usual name for that scale: the
turning line. An empty tick list draws it bare, as above.

Solving propagates: any stage with exactly one unknown is read, its answer becomes known,
and the pass repeats. So inputs can be given in whatever combination closes the chain, and
a combination that cannot is refused with the variables still missing.

The one real constraint is that a shared scale must be the *same* scale in both stages.
Both stages above are multiplicative, so both place `t` logarithmically. Splitting
`V = (L + W) × H` the same way would not work — one stage would want `t` linear and the
other logarithmic — and `compound` says so rather than drawing a chart that lies.

### When a pivot will not do: nets

A pivot works only when two variables collapse into one intermediate. When they do not —
the general quadratic `a·x² + b·x + c = 0` is the classic case — one of the three points the
line crosses has to carry two variables at once. That is a **net**: a mesh of curves,
`(u, v) -> (x, y)`.

```julia
# a·x² + b·x = -c gives α = x², β = x, γ = -c, hence the surface below.
nom = net_nomogram(
    Net(Variable(:x; domain = (0.5, 5)), Variable(:c; domain = (0, 4)),
        (x, c) -> (1 / (1 + x), -c / (x * (1 + x)))),
    Variable(:a; domain = (0.5, 3)) => a -> (0.0, a),
    Variable(:b; domain = (-6, 0)) => b -> (1.0, b);
    relation = (x, c, a, b) -> a * x^2 + b * x + c,
)
[sol.x for sol in evaluate_all(nom; a = 1, b = -3, c = 2)]   # [1.0, 2.0]
```

The reader lays a line from `a` to `b`, finds where it crosses the strand for their `c`, and
reads `x` off the strand through that point. Any of the four variables may be the unknown.
The mesh is drawn lighter than the scales, and each family of strands is labelled at the
end where it fans out rather than where it crowds together.

## Reshaping, for free

A projective map sends straight lines to straight lines, so it preserves the one thing a
nomogram asserts. `project` reshapes a chart **without changing a single answer it gives**:

```julia
corners = [(0, 0), (70, 0), (70, 100), (0, 100)]
skew = project(nom, corners => [(14, 0), (56, 0), (70, 100), (0, 100)])

evaluate(nom;  R1 = 4, R2 = 6).R   # 2.4
evaluate(skew; R1 = 4, R2 = 6).R   # 2.4 — the same chart, differently drawn
```

Four point correspondences fix the map, so you say where the corners should go rather than
writing out nine coefficients. That makes reshaping a design tool rather than a risk: a
chart cramped at one end can be pulled open, and nothing needs re-deriving.

## Examples

- [`examples/parallel_resistance.jl`](examples/parallel_resistance.jl) — resistors in
  parallel, the classic reciprocal chart.
- [`examples/compound_growth.jl`](examples/compound_growth.jl) — compound growth, with a
  custom transform and hand-picked ticks for an axis that compresses hard.
- [`examples/speed_distance_time.jl`](examples/speed_distance_time.jl) — an N-chart for a
  ratio, on linear scales.
- [`examples/quadratic.jl`](examples/quadratic.jl) — roots of `z² + pz + q`, with the
  curved scale derived from the determinant and both roots read off one line.
- [`examples/projective.jl`](examples/projective.jl) — reshaping a chart without changing
  what it computes.
- [`examples/box_volume.jl`](examples/box_volume.jl) — four variables over two stages, with
  a blank turning line.
- [`examples/general_quadratic.jl`](examples/general_quadratic.jl) — four variables that
  no pivot can link, on a net.
- [`examples/from_formula.jl`](examples/from_formula.jl) — compound growth, parallel
  resistors and BMI, each from the formula as written.

## Status

Version 0.1.0, not yet registered. Implemented and tested: building a chart directly from a
formula; the parallel-axis, N-chart and general determinant constructions; compound and net
charts for four or more variables; projective reshaping; the solver; both renderers; and
the tick engine.

Not yet implemented: taking an equation in implicit form, such as `1/R ~ 1/R₁ + 1/R₂`,
rather than a function already solved for its result. That is where a Symbolics.jl front
end would come in, alongside closed-form transforms where they exist and automatic
derivation of curved and net layouts for formulas that are not separable.
