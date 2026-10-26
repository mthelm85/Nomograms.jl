```@meta
CurrentModule = Nomograms
```

# Gallery

This page first builds a chart by hand from plain curves, then shows a range of other
charts, each with the code that draws it.

## A chart built piece by piece

Every chart in this package is a set of [`Scale`](@ref)s: curves that take a value to a
point on the page. The constructions, such as [`additive_nomogram`](@ref), work out those
curves for you. You can also write them yourself, and the result solves, checks and draws
like any other chart.

Here is a classic chart that none of the constructions produces. Take three straight
scales radiating from one point at 30°, 90° and 150°, each graduated evenly with the same
unit. A line through `a` on one outer scale and `b` on the other crosses the middle scale
at `c`, where

```math
\frac{1}{a} + \frac{1}{b} = \frac{1}{c}
```

This is the formula for resistors in parallel, or the thin-lens formula, read from evenly
spaced scales. (For points at distances ``r_i`` along lines through one point at angles
``\theta_i``, collinearity is ``\sum \sin(\theta_j - \theta_k)/r_i = 0`` over the three
cyclic permutations; at these angles that is the relation above.)

**1. The variables.** Each input needs a domain. The result's range follows from the
geometry, and here ``c`` runs from 0.5 to 5:

```@example star
using Nomograms

va = Variable(:a; domain = (1, 10), label = "R₁ (Ω)")
vb = Variable(:b; domain = (1, 10), label = "R₂ (Ω)")
vc = Variable(:c; domain = (0.5, 5), label = "R (Ω)")
nothing # hide
```

**2. The scales.** A [`Scale`](@ref) is a variable, its domain, a curve `value -> (x, y)`,
and the side its labels go on. Each ray here has ten chart units per ohm, and the centre is
placed so that the chart sits at the left of its canvas:

```@example star
centre = (90.0, 0.0)
ray(θ) = r -> (centre[1] + 10r * cosd(θ), centre[2] + 10r * sind(θ))

scales = [
    Scale(va, (1, 10), ray(30), :right),    # labels below the right-hand ray
    Scale(vb, (1, 10), ray(150), :left),    # and below the left-hand one
    Scale(vc, (0.5, 5), ray(90), :right),
]
nothing # hide
```

**3. The chart.** Give the scales a canvas and, optionally, the relation they encode, as a
function that is zero on every solution:

```@example star
star = Nomogram(
    scales;
    width = 180,
    height = 100,
    title = "Resistors in parallel",
    relation = (a, b, c) -> 1 / a + 1 / b - 1 / c,
)
nothing # hide
```

**4. Check it.** The package derived none of this, so check that the geometry encodes the
relation. [`max_relation_residual`](@ref) solves the chart over a grid of inputs and
measures each answer against `relation`:

```@example star
max_relation_residual(star)
```

```@example star
evaluate(star; a = 4, b = 6)
```

**5. Draw it**, with a reading:

```@example star
drawing(star; isopleths = [isopleth(star; a = 4, b = 6, label = "4 ∥ 6 = 2.4")])
```

**6. Style it.** A [`NomogramStyle`](@ref) controls how everything is drawn. Here the tick
labels are turned to run along their scales:

```@example star
style = NomogramStyle(; ticklabelrotation = :tangent)
drawing(star; style, isopleths = [isopleth(star; a = 4, b = 6)])
```

To save it, use `save_svg("star.svg", star; style)`, or `nomogramplot(star; style)` for a
Makie figure. The same steps build any chart you can describe as curves. The
[determinant method](20-constructions.md) is a systematic way to find the curves.

## Logarithmic scales

The most common nomogram is a product, drawn on three parallel logarithmic scales by
[`multiplicative_nomogram`](@ref). This one is Ohm's law, with current from 1 mA to 10 A:

```@example gallery
using Nomograms

ohm = multiplicative_nomogram(
    Variable(:I; domain = (0.001, 10), label = "current (A)"),
    Variable(:R; domain = (1, 1000), label = "resistance (Ω)"),
    Variable(:V; label = "voltage (V)");
    title = "Ohm's law",
)
drawing(ohm; isopleths = [isopleth(ohm; I = 0.25, R = 47, label = "11.75 V")])
```

## An N-chart

[`ratio_nomogram`](@ref) joins two parallel scales with a diagonal. The outer scales stay
linear, and here the labels on the diagonal are turned to run along it:

```@example gallery
pace = ratio_nomogram(
    Variable(:distance; domain = (1, 42), label = "distance (km)"),
    Variable(:time; domain = (5, 240), label = "time (min)"),
    Variable(:speed; label = "speed (km/min)");
    title = "Pace",
)
drawing(
    pace;
    style = NomogramStyle(; ticklabelrotation = :tangent),
    isopleths = [isopleth(pace; distance = 10, time = 50, label = "10 km in 50 min")],
)
```

## A curved scale

Some relations need a curved scale. The roots of ``z^2 + pz + q = 0`` lie on one, and a
line through ``p`` and ``q`` crosses it once for each positive root:

```@example gallery
roots = determinant_nomogram(
    Variable(:p; domain = (-6, 0), label = "p") => p -> (0.0, p),
    Variable(:q; domain = (0, 6), label = "q") => q -> (1.0, q),
    Variable(:z; domain = (0.2, 5), label = "root z") => z -> (1 / (1 + z), -z^2 / (1 + z));
    relation = (p, q, z) -> z^2 + p * z + q,
    title = "z² + pz + q = 0",
)
drawing(roots; isopleths = [isopleth(roots; p = -3, q = 2, label = "roots 1 and 2")])
```

## A net

When two variables have to share one point, that point lies on a [`Net`](@ref): two
families of curves, each named past the end of its row of values.

```@example gallery
timber = net_nomogram(
    Net(
        Variable(:x; domain = (0.5, 5.0), label = "thickness ratio (x)"),
        Variable(:c; domain = (0.0, 4.0), label = "deflection limit (c)"),
        (x, c) -> (1 / (1 + x), -c / (x * (1 + x))),
    ),
    Variable(:a; domain = (0.5, 3.0), label = "load factor (a)") => a -> (0.0, a),
    Variable(:b; domain = (-6.0, 0.0), label = "span factor (b)") => b -> (1.0, b);
    relation = (x, c, a, b) -> a * x^2 + b * x + c,
    title = "ax² + bx + c = 0",
)
drawing(timber; isopleths = [isopleth(timber; a = 1, b = -3, c = 2)])
```

## A chain

Three inputs need two alignment lines with an unmarked turning line between them.
[`auto_nomogram`](@ref) finds the chain from the formula. This example uses units, so the
pressure comes out in kPa whatever units the inputs are given in:

```@example gallery
using Unitful

R = 8.314u"J/(mol*K)"
gas = auto_nomogram(
    (n, T, V) -> n * R * T / V,
    Variable(:n; domain = (0.5u"mol", 5u"mol"), label = "amount"),
    Variable(:T; domain = (250u"K", 400u"K"), label = "temperature"),
    Variable(:V; domain = (10u"L", 100u"L"), label = "volume"),
    Variable(:P; unit = u"kPa", label = "pressure");
    title = "Ideal gas",
)
drawing(gas; isopleths = isopleth(gas; n = 1u"mol", T = 300u"K", V = 24u"L", label = "104 kPa"))
```

## Converging scales

A projective map keeps straight lines straight, so it can reshape a chart without changing
its answers. [`project`](@ref) takes the corners of the canvas and where they should go:

```@example gallery
parallel = additive_nomogram(
    Variable(:R1, ReciprocalScale(); domain = (1, 10), label = "R₁ (Ω)"),
    Variable(:R2, ReciprocalScale(); domain = (1, 10), label = "R₂ (Ω)"),
    Variable(:R, ReciprocalScale(); label = "R (Ω)");
    title = "Converging",
)
corners = [(0, 0), (70, 0), (70, 100), (0, 100)]
converging = project(parallel, corners => [(0, 0), (70, 0), (50, 100), (20, 100)])
drawing(converging; isopleths = [isopleth(converging; R1 = 4, R2 = 6, label = "2.4")])
```

## A fan

Take the same idea further, until the scales meet at a point beyond the chart, and they
become spokes. [`fan`](@ref) does this for any chart, including compound ones:

```@example gallery
fanned = fan(ohm; spread = 80)
drawing(fanned; isopleths = [isopleth(fanned; I = 0.25, R = 47, label = "11.75 V")])
```

## A ruler

Fraction ticks are graded by importance, as on a ruler: wholes are longest, then halves,
quarters, eighths and sixteenths. Writing the domains as `Rational`s turns them on:

```@example gallery
bolt = additive_nomogram(
    Variable(:grip; domain = (1//4, 3), label = "grip (in)"),
    Variable(:nut; domain = (1//8, 3//4), maxden = 16, label = "nut + washers (in)"),
    Variable(:length; tickstyle = :fraction, label = "bolt length (in)");
    title = "Bolt length",
)
drawing(bolt; isopleths = [isopleth(bolt; grip = 1 + 1//4, nut = 3//8, label = "1 5/8 in")])
```

## A Wulff net

Crystallographers and structural geologists plot directions in space on a Wulff net: the
stereographic projection of a sphere's latitude and longitude lines, with the poles at the
top and bottom. It is a [concurrency chart](35-concurrency-charts.md) on a plane of its
own, like the Smith chart, built from the same three keywords: `axes = false` for the
plane, `rings` for the azimuth scale around the edge, and `labelat` to label latitude down
the central meridian and longitude along the equator.

```@example gallery
# The point of the unit sphere that projects to (u, v).
function sphere(u, v)
    r2 = u^2 + v^2
    return (2u / (1 + r2), 2v / (1 + r2), (1 - r2) / (1 + r2))
end
latitude(u, v) = asind(clamp(sphere(u, v)[2], -1, 1))
longitude(u, v) = (s = sphere(u, v); atand(s[1], s[3]))

wulff = concurrency_chart(
    :u => (-1, 1),
    :v => (-1, 1),
    Variable(:lat; domain = (-90, 90), ticks = -80:10:80, minorticks = -88:2:88,
        label = "latitude (°)") => latitude,
    Variable(:lon; domain = (-90, 90), ticks = -80:10:80, minorticks = -88:2:88,
        label = "longitude (°)") => longitude;
    axes = false,
    region = (u, v) -> u^2 + v^2 <= 1,
    rings = [
        Variable(:azimuth; domain = (0, 360), ticks = 0:30:330, minorticks = 0:10:360,
            label = "azimuth (°)") => a -> (1.06sind(a), 1.06cosd(a)),
    ],
    labelat = (lat = :lon => 0, lon = :lat => 0),
    title = "Wulff net",
)
drawing(wulff; crossings = [crossing(wulff; lat = 30, lon = 40, label = "30°, 40°")])
```

The grid is every 2°, but the 2° meridians stop short of the poles, where they converge,
and only every tenth one carries on. Printed nets are drawn the same way. Here it comes
from the spacing rule that thins every concurrency chart's curves, with nothing specific
to this chart.

## Colours

[`svg`](@ref) accepts colours for the background, axes, text and reading lines, so a
chart can match the colour scheme of the page it appears on:

```@example gallery
struct SVGText; s::String; end # hide
Base.show(io::IO, ::MIME"image/svg+xml", x::SVGText) = print(io, x.s) # hide
dark = svg(
    ohm;
    isopleths = [isopleth(ohm; I = 0.25, R = 47, label = "11.75 V")],
    background = "#14213d",
    axiscolor = "#e5e5e5",
    textcolor = "#e5e5e5",
    isoplethcolors = ("#fca311",),
)
SVGText(dark) # hide
```
