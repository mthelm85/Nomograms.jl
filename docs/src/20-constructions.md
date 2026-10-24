```@meta
CurrentModule = Nomograms
```

# Building charts

A nomogram asserts one thing: that a point taken from each scale lies on a single straight
line. So every chart in this package is a set of parametric curves `value -> (x, y)`
([`Scale`](@ref)), and each construction below is a way of choosing curves that makes that
collinearity agree with a relation. [How it works](60-how-it-works.md) follows a chart
from construction to picture in detail.

## Parallel scales

A three-variable relation fits on three parallel scales exactly when it can be separated
as

```math
f_1(u) + f_2(v) = f_3(w)
```

Each `fᵢ` is a [`ScaleTransform`](@ref): [`LinearScale`](@ref), [`LogScale`](@ref),
[`ReciprocalScale`](@ref), [`PowerScale`](@ref), [`ComposedScale`](@ref), or
[`FunctionScale`](@ref) for anything else. A transform carries its own inverse, which is
what lets the result axis be derived exactly instead of searched for.
[`additive_nomogram`](@ref) builds the chart; see its docstring for the geometry.

Where the middle axis goes is not a free choice: it falls between the outer two, nearer the
one whose transform spans more. When one range is far wider than the other, the middle axis
would sit almost on top of that neighbour and its labels would read as the neighbour's. So
the narrower axis is shortened and centred until the middle axis is at least a fifth of the
way across. Readings are unchanged, since the middle axis moves with it. Set
`max_modulus_ratio` to trade the narrow axis's length against that spacing, or `Inf` for
full-height axes whatever their ranges.

Product relations separate under a logarithm, which [`multiplicative_nomogram`](@ref)
applies for you:

```@example build
using Nomograms

cost = multiplicative_nomogram(
    Variable(:rate; domain = (1, 100), label = "rate (\$/h)"),
    Variable(:hours; domain = (1, 100), label = "hours"),
    Variable(:cost; label = "cost (\$)");
    title = "Labour cost",
)
drawing(cost; isopleths = [isopleth(cost; rate = 7.5, hours = 12, label = "\$90")])
```

## Ratios, and the N-chart

Not every relation wants parallel scales. [`ratio_nomogram`](@ref) builds an N-chart (also
called a Z-chart) for `f₁(u) / f₂(v) = f₃(w)`: two vertical scales joined by a diagonal,
with the result read off the diagonal. Doerfler's "The Lost Art of Nomography" derives the
classic construction (see the [Bibliography](96-references.md)).

```@example build
speed = ratio_nomogram(
    Variable(:distance; domain = (1, 10), label = "distance (km)"),
    Variable(:time; domain = (1, 10), label = "time (h)"),
    Variable(:speed; label = "speed (km/h)");
    title = "Speed",
)
drawing(speed; isopleths = [isopleth(speed; distance = 6, time = 3, label = "2 km/h")])
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

## Anything else, from the determinant

The two named constructions are particular solutions of the collinearity condition.
[`determinant_nomogram`](@ref) takes any other one you can write down: give each variable a
curve, and the combined bounding box is measured and mapped onto the canvas for you. This
is the determinant method of the nomography textbooks; Doerfler's overview is a good
worked introduction to it.

Deriving a curve is mechanical once the first two scales are fixed. Put `p` on `x = 0` at
height `p` and `q` on `x = 1` at height `q`; a line between them crosses `x = x₃` at
`p + x₃(q - p)`, so a third scale is met exactly when `(1 - x₃)·p + x₃·q = y₃`. To make
that say `z·p + q + z² = 0`, match coefficients:

```math
\frac{1 - x_3}{z} = \frac{x_3}{1} = \frac{-y_3}{z^2}
```

giving `x₃ = 1/(1 + z)` and `y₃ = -z²/(1 + z)` — a curved scale for the roots of a
quadratic.

```@example build
quad = determinant_nomogram(
    Variable(:p; domain = (-6, 0)) => p -> (0.0, p),
    Variable(:q; domain = (0, 6)) => q -> (1.0, q),
    Variable(:z; domain = (0.2, 5), label = "root z") => z -> (1 / (1 + z), -z^2 / (1 + z));
    relation = (p, q, z) -> z^2 + p * z + q,
    title = "z² + pz + q = 0",
)
roots = evaluate_all(quad; p = -3, q = 2)
drawing(quad; isopleths = [isopleth(quad, roots[1]; label = "z² − 3z + 2 = 0")])
```

```@example build
[sol.z for sol in roots]
```

Two answers, because a quadratic has two roots and the chart shows why: one alignment line,
meeting the curved scale twice. This curved-scale chart for the quadratic is a classic;
Doerfler draws it after Adams (1964) and Otto (1963).

Collinearity survives any affine map, so fitting the curves to the canvas cannot change
what the chart computes — including the anisotropic fit used by default to fill it. Pass
`preserve_aspect = true` to scale both axes alike and centre the result.

Getting a determinant form right is the hard part of designing a chart by hand, and a chart
that is subtly wrong still draws perfectly well. [`max_relation_residual`](@ref) solves on
the geometry over a grid of inputs and checks each answer against the algebra:

```@example build
max_relation_residual(quad)
```

It returns `NaN` rather than a flattering zero when there is no relation to check, or when
no input in the grid produced a solution.

## Reshaping, for free

An affine map preserves collinearity, which is why fitting curves to the canvas is safe.
*Projective* maps preserve it too, and they are a much larger family — they can turn
parallel scales into converging ones and redistribute spacing along an axis. Nomographers
have long reshaped charts this way; Doerfler's overview closes with these transformations.

So [`project`](@ref) reshapes a chart without changing a single answer it gives:

```@example build
parallel = additive_nomogram(
    Variable(:R1, ReciprocalScale(); domain = (1, 10), label = "R₁ (Ω)"),
    Variable(:R2, ReciprocalScale(); domain = (1, 10), label = "R₂ (Ω)"),
    Variable(:R, ReciprocalScale(); label = "R (Ω)"),
)
corners = [(0, 0), (70, 0), (70, 100), (0, 100)]
skew = project(parallel, corners => [(14, 0), (56, 0), (70, 100), (0, 100)])
drawing(skew; isopleths = [isopleth(skew; R1 = 4, R2 = 6, label = "2.4")])
```

```@example build
(evaluate(parallel; R1 = 4, R2 = 6).R, evaluate(skew; R1 = 4, R2 = 6).R)
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

## Fans

Parallel lines meet at a point at infinity, and a projective map can bring that point onto
the page. Every axis then runs through it, and the chart becomes a fan of spokes —
[`fan`](@ref) does this:

```@example build
title = "Resistors in parallel"
fanned = fan(additive_nomogram(
    Variable(:R1, ReciprocalScale(); domain = (1, 10), label = "R₁ (Ω)"),
    Variable(:R2, ReciprocalScale(); domain = (1, 10), label = "R₂ (Ω)"),
    Variable(:R, ReciprocalScale(); label = "R (Ω)");
    title,
))
drawing(fanned; isopleths = [isopleth(fanned; R1 = 4, R2 = 6, label = "4 ∥ 6 = 2.4")])
```

It is the same chart, reading the same answers, and works on compound charts and nets too.
Two things follow from the geometry. Straight lines stay straight, so the spokes end on a
straight chord rather than an arc: at the default 90° spread the tips line up level, and
wider fans leave the middle spokes shorter than the outer ones. And the map squeezes whatever lies
near the hub, so `fan` lays the chart out both ways round and keeps whichever labels more
values. Here the reciprocal scales' crowded ends — 4 to 10 — went to the rim.
