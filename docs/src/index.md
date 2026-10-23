```@meta
CurrentModule = Nomograms
```

# Nomograms

Build, solve and render nomograms — alignment charts that turn an equation into a
graphical calculator, where any straight line laid across the scales picks out a solution.

Source: [Nomograms.jl](https://github.com/mthelm85/Nomograms.jl).

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/mthelm85/Nomograms.jl")
```

## A first chart

Two resistors in parallel obey `1/R = 1/R₁ + 1/R₂`, which is already in the additive form
a parallel-axis nomogram needs; the transform is simply `1/x`.

```@example intro
using Nomograms

nom = additive_nomogram(
    Variable(:R1, ReciprocalScale(); domain = (1, 10), label = "R₁ (Ω)"),
    Variable(:R2, ReciprocalScale(); domain = (1, 10), label = "R₂ (Ω)"),
    Variable(:R,  ReciprocalScale(); label = "R (Ω)");
    title = "Resistors in parallel",
)
```

`R` declares no domain. The length and position of the middle axis follow from the outer
two rather than being free parameters, so the construction derives the range `R` can
reach:

```@example intro
nom[:R].domain
```

Solve it in either direction:

```@example intro
evaluate(nom; R1 = 4, R2 = 6)
```

```@example intro
evaluate(nom; R1 = 4, R = 2.4)
```

And render it, either with the built-in dependency-free SVG writer:

```julia
save_svg("parallel.svg", nom; isopleths = [isopleth(nom; R1 = 4, R2 = 6)])
```

or through Makie, once a backend is loaded:

```julia
using CairoMakie
save("parallel.pdf", nomogramplot(nom; isopleths = [isopleth(nom; R1 = 4, R2 = 6)]))
```

## Starting from a formula

The constructions above ask you to rewrite a formula into a separable form and choose its
transforms. [`auto_nomogram`](@ref) does that part, from the function as you would
ordinarily write it:

```@example intro
growth = auto_nomogram((r, n) -> (1 + r / 100)^n, :r => (1, 20), :n => (1, 30), :growth)
evaluate(growth; r = 7, growth = 2)
```

### When a formula has a chart at all

A formula `z = F(x, y)` can be drawn on three straight parallel scales exactly when there
are transforms with ``h(z) = f(x) + g(y)``. Differentiating that with respect to each
argument gives

```math
h'(z)\,F_x = f'(x), \qquad h'(z)\,F_y = g'(y)
\quad\Longrightarrow\quad
\frac{F_x}{F_y} = \frac{f'(x)}{g'(y)}
```

so ``\ln|F_x / F_y|`` must split into a function of ``x`` plus a function of ``y``. The
converse holds too: if it splits, take ``f'`` and ``g'`` from the two parts; then
``\nabla(f + g)`` is everywhere parallel to ``\nabla F``, so ``f + g`` is constant exactly
where ``F`` is, which makes it a function of ``F`` — and that function is ``h``. The test
is classical in nomography, but the argument above is all it rests on here.

[`separability_defect`](@ref) applies it numerically. Rather than take a mixed second
derivative of a logarithm, which is noisy, it checks the equivalent identity

```math
L(x, y) - L(x, y_0) - L(x_0, y) + L(x_0, y_0) = 0, \qquad L = \ln|F_x / F_y|
```

which needs only first derivatives:

```@example intro
separability_defect((x, y) -> x + y + x * y, (0, 5), (0, 5)),   # (1+x)(1+y) - 1
separability_defect((x, y) -> x^2 + x * y + y, (0.5, 5), (0.5, 5))
```

### Recovering the transforms

When the test passes, the same identity hands the transforms back. Fixing ``y`` at a
reference value ``y_0`` gives ``f'(x)`` up to a constant, fixing ``x`` at ``x_0`` gives
``g'(y)``, and numerical integration turns each into a [`TabulatedScale`](@ref). ``h`` is
tabulated by walking from the corner of the domain that minimises ``F`` to the one that
maximises it, along which ``f`` and ``g`` change the same way, so ``F`` never has to be
inverted.

Each axis then gets linear or decade-based ticks according to how its transform spreads its
values: if the geometric mean of the domain lands nearer the middle of the axis than the
arithmetic mean, the axis is logarithmic in character and gets decade ticks.

### Nothing on trust

The finished chart is *solved* over a grid of inputs and every answer compared with `F`
itself. A chart that disagrees by more than `rtol` of the result's range is not returned —
the call throws instead, quoting the disagreement it measured. The chart's `relation` is
your function, so [`max_relation_residual`](@ref) measures it against what you wrote rather
than anything derived along the way.

A formula with no straight-scale chart gets an error saying so, pointing at
[`determinant_nomogram`](@ref) and [`net_nomogram`](@ref), which can still handle it given a
form derived by hand. So does a formula that is not monotone in some variable over its
domain, since a straight scale cannot fold back on itself.

## How it works

A nomogram asserts one thing: that a point taken from each scale lies on a single straight
line. Written as a determinant, that is d'Ocagne's criterion,

```math
\begin{vmatrix}
x_1 & y_1 & 1\\
x_2 & y_2 & 1\\
x_3 & y_3 & 1
\end{vmatrix} = 0
```

so this package stores every scale as a parametric curve `value -> (x, y)` ([`Scale`](@ref))
and does all of its arithmetic by asking where a line crosses a curve
([`align_residual`](@ref), [`evaluate`](@ref)). One consequence is that forward and inverse
problems are the same operation — any variable may be the unknown. Another is that a
layout the built-in constructions do not cover needs no new solver: hand-build the scales
and everything downstream still works.

### More than one crossing

A straight line can meet a *curved* scale more than once, and each crossing is a genuine
solution of the relation. The usual shortcut — checking the sign of the residual at the two
ends of the domain and bisecting — gets this wrong twice over: it reports "no solution" for
an even number of crossings, because both ends lie on the same side of the line, and it
returns an arbitrary one of an odd number without saying so.

So [`solve_all_on`](@ref) searches the whole domain. [`evaluate`](@ref) still returns a
single answer, but throws rather than choosing when there are several, naming the values it
could not choose between; [`evaluate_all`](@ref) returns every solution:

```julia
sols = evaluate_all(nom; u = 20, v = 20)
isos = [isopleth(nom, sol) for sol in sols]
```

A straight scale makes the residual an affine function of position, hence monotone, so
there is always at most one root: parallel-layout nomograms never encounter this. The cost
of searching rather than bracketing is around 120 µs per `evaluate` instead of a few, which
is immaterial for reading charts.

### When is a relation chartable?

A three-variable relation fits on three parallel scales exactly when it can be separated
as

```math
f_1(u) + f_2(v) = f_3(w)
```

Each `fᵢ` is a [`ScaleTransform`](@ref): [`LinearScale`](@ref), [`LogScale`](@ref),
[`ReciprocalScale`](@ref), [`PowerScale`](@ref), [`ComposedScale`](@ref), or
[`FunctionScale`](@ref) for anything else. A transform carries its own inverse, which is
what lets the result axis be derived exactly instead of searched for.

Product relations separate under a logarithm, which
[`multiplicative_nomogram`](@ref) applies for you:

```@example intro
prod = multiplicative_nomogram(
    Variable(:rate;  domain = (1, 100)),
    Variable(:hours; domain = (1, 100)),
    Variable(:cost),
)
evaluate(prod; rate = 7.5, hours = 12)
```

See [`additive_nomogram`](@ref) for the geometry of the construction, including why the
middle axis spans exactly the same height as the outer two.

### Ratios, and the N-chart

Not every relation wants parallel scales. [`ratio_nomogram`](@ref) builds an N-chart (also
called a Z-chart) for `f₁(u) / f₂(v) = f₃(w)`: two vertical scales joined by a diagonal,
with the result read off the diagonal.

```@example intro
speed = ratio_nomogram(
    Variable(:distance; domain = (1, 10)),
    Variable(:time; domain = (1, 10)),
    Variable(:speed),
)
evaluate(speed; distance = 6, time = 3)
```

`multiplicative_nomogram` charts ratios as well, by way of logarithms. The N-chart earns
its place by needing none: the outer scales stay *linear*, so precision is even along them
and a reader can interpolate between ticks by eye.

Two details are worth knowing. The textbook derivation gives both outer scales the same
scale factor so that the ratio falls out directly; here the diagonal's own scale absorbs
the ratio of the two factors instead, so each outer axis is fitted to the canvas on its
own. And the diagonal is a *projective* scale, `s = k·r / (k·r + 1)`, so a ratio spanning
decades is spread across it roughly logarithmically — which is why the result axis takes
decade-based tick values in that case. Picking them linearly would offer nothing between
0.1 and 2, most of the diagonal. Pass `tickstyle = :linear` on the [`Variable`](@ref) to
override.

### Anything else, from the determinant

The two named constructions are particular solutions of the collinearity condition.
[`determinant_nomogram`](@ref) takes any other one you can write down: give each variable a
curve, and the combined bounding box is measured and mapped onto the canvas for you.

Deriving a curve is mechanical once the first two scales are fixed. Put `p` on `x = 0` at
height `p` and `q` on `x = 1` at height `q`; a line between them crosses `x = x₃` at
`p + x₃(q - p)`, so a third scale is met exactly when `(1 - x₃)·p + x₃·q = y₃`. To make
that say `z·p + q + z² = 0`, match coefficients:

```math
\frac{1 - x_3}{z} = \frac{x_3}{1} = \frac{-y_3}{z^2}
```

giving `x₃ = 1/(1 + z)` and `y₃ = -z²/(1 + z)` — a curved scale for the roots of a
quadratic.

```@example intro
quad = determinant_nomogram(
    Variable(:p; domain = (-6, 0)) => p -> (0.0, p),
    Variable(:q; domain = (0, 6)) => q -> (1.0, q),
    Variable(:z; domain = (0.2, 5)) => z -> (1 / (1 + z), -z^2 / (1 + z));
    relation = (p, q, z) -> z^2 + p * z + q,
)
[sol.z for sol in evaluate_all(quad; p = -3, q = 2)]
```

Two answers, because a quadratic has two roots and the chart shows why: one alignment line,
meeting the curved scale twice.

Collinearity survives any affine map, so fitting the curves to the canvas cannot change
what the chart computes — including the anisotropic fit used by default to fill it. Pass
`preserve_aspect = true` to scale both axes alike and centre the result.

Getting a determinant form right is the hard part of designing a chart by hand, and a chart
that is subtly wrong still draws perfectly well. [`max_relation_residual`](@ref) solves on
the geometry over a grid of inputs and checks each answer against the algebra:

```@example intro
max_relation_residual(quad)
```

It returns `NaN` rather than a flattering zero when there is no relation to check, or when
no input in the grid produced a solution.

### More than three variables

One alignment line can relate only three scales. A fourth variable means two stages sharing
an intermediate, and [`compound`](@ref) positions them so the shared scale coincides.

```@example intro
stage1 = multiplicative_nomogram(
    Variable(:L; domain = (1, 10)),
    Variable(:W; domain = (1, 10)),
    Variable(:t; label = "", ticks = Float64[]),
)
stage2 = multiplicative_nomogram(
    Variable(:t; domain = stage1[:t].domain, label = "", ticks = Float64[]),
    Variable(:H; domain = (1, 10)),
    Variable(:V);
    width = 130,
)
box = compound(stage1, stage2; title = "V = L × W × H")
evaluate(box; L = 3, W = 4, H = 2)
```

The reader lays a line across the first stage, notes where it crosses the shared scale,
then pivots about that point and reads the second — which is why that scale is called the
turning line. Giving it an empty tick list draws it bare, as here; leave the ticks on and
it carries the intermediate value, which is worth doing when the intermediate means
something.

Solving propagates rather than following a fixed direction: any stage with exactly one
unknown is read, its answer becomes known, and the pass repeats. So inputs may be given in
whatever combination closes the chain,

```@example intro
evaluate(box; L = 3, W = 4, V = 24).H
```

and a combination that cannot close it is refused, naming what stayed unknown. Supplying
`t`, `H` and `V`, for instance, pins the product `L × W` but neither factor.

Alignment is a similarity — rotate, scale, translate — carrying the later stage's shared
scale onto the placed one, checked along the whole domain rather than only at the two
points that fixed it. Being affine, it cannot change what either stage computes.

!!! note "The shared scale must be the same scale in both stages"
    Both stages above are multiplicative, so both place `t` logarithmically and the two
    coincide. Splitting `V = (L + W) × H` the same way would not work: the first stage
    puts `t` on a linear scale and the second on a logarithmic one, and no amount of
    moving the stages reconciles that. `compound` checks and says so rather than drawing
    a chart that lies.

Stages generally interleave rather than sitting side by side, because a parallel-axis
construction puts its result *between* its inputs, so both stages straddle the pivot.
Widening the later stage's canvas spreads it clear, and [`project`](@ref) will reshape the
result if the arrangement is still awkward.

### When a pivot will not do: nets

A pivot works only when the first two variables collapse into a single intermediate. When
they do not, one of the three points an alignment line crosses must carry two variables at
once. That is a [`Net`](@ref): a two-parameter surface `(u, v) -> (x, y)`, drawn as a mesh.

Deriving one follows the same pattern as a curved scale. With `a` on `x = 0` at height
`f(a)` and `b` on `x = 1` at height `g(b)`, a net point is crossed exactly when

```math
\alpha f(a) + \beta g(b) = \gamma, \qquad
x_3 = \frac{\beta}{\alpha + \beta}, \quad y_3 = \frac{\gamma}{\alpha + \beta}
```

so any relation whose two net variables enter as the coefficients has a net. The general
quadratic `a·x² + b·x + c = 0`, written `a·x² + b·x = -c`, has `α = x²`, `β = x`, `γ = -c`:

```@example intro
gq = net_nomogram(
    Net(
        Variable(:x; domain = (0.5, 5)),
        Variable(:c; domain = (0, 4)),
        (x, c) -> (1 / (1 + x), -c / (x * (1 + x))),
    ),
    Variable(:a; domain = (0.5, 3)) => a -> (0.0, a),
    Variable(:b; domain = (-6, 0)) => b -> (1.0, b);
    relation = (x, c, a, b) -> a * x^2 + b * x + c,
)
[sol.x for sol in evaluate_all(gq; a = 1, b = -3, c = 2)]
```

Which problem [`evaluate`](@ref) solves depends on what is unknown. With both net variables
known the net contributes a single point and it is the ordinary three-point alignment. With
one unknown, the line is already fixed by the two scales, and the answer is where it
crosses the strand for the known net variable — a curve, which the line may cross twice.
Here it does, because a quadratic has two roots.

The mesh is genuinely two-dimensional only when `(α : β : γ)` sweeps a region. When it
sweeps just a curve the net collapses onto a line, and the relation wanted a pivot instead —
which is a useful check on whether you have reached for the right construction.

Strands are scales of `kind = :grid`, drawn lighter than the scales a reader takes values
from. Each family is labelled at whichever end it fans out, since at the end where strands
converge the labels pile up and get nudged next to the wrong curve; set `titleend` on a
[`Variable`](@ref) to choose by hand.

### Reshaping, for free

An affine map preserves collinearity, which is why fitting curves to the canvas is safe.
*Projective* maps preserve it too, and they are a much larger family — they can turn
parallel scales into converging ones and redistribute spacing along an axis.

So [`project`](@ref) reshapes a chart without changing a single answer it gives:

```@example intro
corners = [(0, 0), (70, 0), (70, 100), (0, 100)]
skew = project(nom, corners => [(14, 0), (56, 0), (70, 100), (0, 100)])
(evaluate(nom; R1 = 4, R2 = 6).R, evaluate(skew; R1 = 4, R2 = 6).R)
```

[`projective_transform`](@ref) builds the map from four point correspondences, which is
how you want to say it in practice: name the corners of a quadrilateral and where they
should end up, rather than writing out nine coefficients. A map whose bottom row is
`[0 0 1]` is affine — the familiar translate, scale and shear — and everything above that
is genuine perspective.

The result is re-fitted to the canvas, since a projective map takes the bounding box
somewhere arbitrary; composing that affine fit with the projective map leaves it
projective, so the guarantee is unaffected.

This is what makes reshaping a design tool rather than a risk. A region where the ticks
pile up can be pulled open and nothing needs re-deriving — and
[`max_relation_residual`](@ref) will confirm it on the reshaped chart.

The one thing to watch for is the map's vanishing line, where the denominator goes to
zero. If it crosses a scale's domain, points there go to infinity; that is reported
against the scale it affects rather than producing a chart full of `Inf`.

### Ticks

Tick values are chosen to be round numbers *in the variable*, then thinned by distance *on
the page* ([`thin_ticks`](@ref)). That split is what keeps a non-uniform axis readable: a
reciprocal or logarithmic scale keeps its labels where it is stretched and drops them
where it is compressed, and the numbers that survive are still round ones. Axis endpoints
and decade boundaries are placed first, so a crowded axis gives up `9` before it gives up
`10`.

Where an axis compresses so hard that automatic selection has nothing to work with, pass
`ticks = [...]` to the [`Variable`](@ref); anything that will not fit is still dropped.

### Placing the labels

Thinning compares a scale only against itself, which is all a parallel layout needs. Bring
scales together — a diagonal, a curve, two axes close on the page — and labels from
*different* scales start to overlap, which no amount of per-axis thinning can see.

So each label carries its anchor point, alignment and rotation, and
[`label_box`](@ref Nomograms.label_box) turns those into the box it actually occupies. The
layout then makes one pass over the whole chart: axis titles are reserved first, then major
tick labels, then minor ones, and a label that cannot fit is dropped while its tick mark
stays. Isopleth labels are nudged along the perpendicular until they are clear instead,
since the caller asked for those by name.

Collisions are not only label against label. An alignment line drawn straight through a
number makes it just as unreadable, so [`segment_hits_box`](@ref Nomograms.segment_hits_box)
tests each label against the lines too and a struck label is dropped, tick mark intact.
Only the span between the outermost crossings counts: the overhang past them is a flourish,
and it falls exactly where the outer scales' labels sit -- including the label for the value
being read. Set `avoid_isopleths = false` to keep every label regardless.

Two related details fall out of the same geometry:

- [`align_for`](@ref Nomograms.align_for) uses *both* components of the direction a label
  is pushed in. While every axis was vertical only the horizontal one mattered; on a curve,
  a label pushed straight down from its tick has to hang from its top edge rather than sit
  centred on the point.
- An axis title is placed off the **end** of its scale, continuing along the tangent. For a
  vertical axis that is straight up from the top, exactly as before; for a diagonal or a
  curve it follows the scale instead of jumping to whichever point is highest.

Set `ticklabelrotation = :tangent` on the [`NomogramStyle`](@ref) to turn tick labels to run
along their scale, which reads better on a steep or curved axis; the default leaves them
upright. `resolve_collisions = false` turns the whole pass off.

### Rendering

[`drawing`](@ref) turns a nomogram into a [`NomogramDrawing`](@ref): axis polylines, tick
positions, normals and labels, and isopleth segments — pure geometry, with nothing about
how to draw it. Both backends are thin walks over that structure, so a
[`NomogramStyle`](@ref) means the same thing in either, and tests can assert on tick
positions without rendering anything.

| | dependency | output |
|---|---|---|
| [`svg`](@ref), [`save_svg`](@ref) | none | SVG |
| [`nomogramplot`](@ref) | any Makie backend | `Figure`, hence PDF/SVG/PNG |

### Styling

The defaults follow the proportions of a classic d3 axis: hairline rules, short ticks and
small type. Everything in [`NomogramStyle`](@ref) is measured in chart units, so at the
default `pxperunit = 4` they land on 1px strokes, 6px major ticks, 3px of label padding and
10px labels, and changing the output size rescales the lot together.

Isopleths get three colours — blue, orange, aqua — and that is deliberate rather than
parsimonious. Every isopleth on a chart is visible at the same time, so all three pairs
have to stay apart under colour-vision deficiency, and no fourth hue holds up once every
pair is in play. The three were checked with a validator rather than by eye: worst pair
ΔE 9.2 under deuteranopia and 24.0 in normal vision, against a ≥8 target.

d3's own `category10` is not used, for a reason that is measurable rather than aesthetic:
its orange and green are ΔE 0.7 apart under protanopia, which is to say indistinguishable.

Past three isopleths the dash pattern advances instead of the hue — a fourth line is a
dashed blue, not a second plain blue — so identity never rests on colour alone. Override
either with `isoplethcolors` and `isoplethdashes`, which both backends accept.

## Reference

See the [Reference](@ref reference) page for the full API.
