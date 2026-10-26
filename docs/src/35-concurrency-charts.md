```@meta
CurrentModule = Nomograms
```

# Concurrency charts

Every chart elsewhere in this package is read with a straight line laid across its scales.
A concurrency chart is read at a point instead. Two variables run along the axes, and every
other variable is drawn as a family of curves, one curve per round value. Find where the
curves of two values you know cross, and read every other value at that point.

Charts of this kind are older than alignment charts. d'Ocagne devised the alignment method
as a way to replace the crossing curves of the charts of his day with points on a line.
They are still the right tool in two cases:

- **Any formula has one.** An alignment chart needs the formula to have a special form,
  which is why [`auto_nomogram`](@ref) sometimes refuses. A concurrency chart only needs
  the formula's value at each point.
- **Many variables can share one sheet.** Any two values fix a point, and you can read as
  many other values there as there are families of curves. A psychrometric chart shows
  seven properties of moist air this way.

The price is precision. You read between curves, so a concurrency chart is usually less
precise than an alignment chart of the same size. Where an alignment chart exists, it is
usually the better choice.

## A first chart

[`concurrency_chart`](@ref) takes the two axis variables and then any number of *fields*,
each written `variable => (x, y) -> value`. A field needs only a formula for its value at
each point. Here is the Moody diagram for friction in pipes. The Colebrook equation that
defines it cannot be solved for the friction factor, but it can be solved for the pipe's
relative roughness, and that is all the chart needs:

```@example conc
using Nomograms

roughness(Re, f) = 3.7 * (10^(-1 / (2sqrt(f))) - 2.51 / (Re * sqrt(f)))

moody = concurrency_chart(
    Variable(:Re, LogScale(); domain = (4e3, 1e8), label = "Reynolds number"),
    Variable(:f, LogScale(); domain = (0.008, 0.1), label = "friction factor"),
    Variable(:roughness; domain = (1e-6, 0.05), label = "relative roughness ε/D") =>
        roughness;
    region = (Re, f) -> roughness(Re, f) >= 0,
    title = "Moody diagram",
)
drawing(moody; readings = [reading(moody; Re = 1e5, roughness = 1e-4, label = "f = 0.0185")])
```

The axes take `LogScale` transforms, so they are logarithmic, as on the printed chart. The
roughness spans more than two decades, so its curves are placed at logarithmic values
automatically. `region` limits the chart to where the roughness is zero or more. The edge
of the region is the curve for a perfectly smooth pipe, which is drawn as part of the frame.

## Reading a chart

[`evaluate`](@ref) takes any two variables and returns all of them:

```@example conc
evaluate(moody; Re = 1e5, roughness = 1e-4)
```

The two values can be two axes, an axis and a field, or two fields. If their curves cross
more than once, `evaluate` refuses and [`evaluate_all`](@ref) returns every crossing.

To mark a reading on the chart, build it with [`reading`](@ref) and pass it to
[`drawing`](@ref) as `readings`. The point is drawn with lines dropped to both axes.
`isopleth` doesn't apply here: an isopleth is a line through an alignment chart, and this
chart is read at a point.

## Several families

A psychrometric chart has temperature and humidity ratio on its axes, and curves of
relative humidity and enthalpy across them. Air can't hold more water than saturates it,
so the chart stops at the saturation line, where the relative humidity is 1:

```@example conc
P = 101.325                                          # kPa, sea level
saturation(T) = 0.61094 * exp(17.625T / (T + 243.04))  # kPa, Magnus formula
vapour(W) = W * P / (622 + W)                        # kPa, W in g/kg

psy = concurrency_chart(
    Variable(:T; domain = (0, 45), label = "dry-bulb temperature (°C)"),
    Variable(:W; domain = (0, 30), label = "humidity ratio (g/kg)"),
    Variable(:RH; domain = (0.1, 1), ticks = 0.1:0.1:1.0, label = "relative humidity") =>
        (T, W) -> vapour(W) / saturation(T),
    Variable(:h; label = "enthalpy (kJ/kg)") =>
        (T, W) -> 1.006T + W / 1000 * (2501 + 1.86T);
    region = (T, W) -> vapour(W) <= saturation(T) * (1 + 1e-9),
    title = "Psychrometric chart",
)
drawing(psy; readings = [reading(psy; T = 30, RH = 0.5, label = "30 °C, 50%")])
```

Any two properties fix the state of the air. Given relative humidity and enthalpy, two
families of curves, the chart finds where they cross:

```@example conc
evaluate(psy; RH = 0.5, h = 60)
```

## From a formula

Given a function first, `concurrency_chart(F, x, y, z)` charts `z = F(x, y)`, in the same
way [`auto_nomogram`](@ref) takes its arguments. It works for any `F`. Here is the formula
`auto_nomogram` refuses:

```@example conc
concurrency_chart((x, y) -> sin(x) + exp(y) + x * y, :x => (0, 2), :y => (0, 1), :z)
```

When `F` passes the separability test that `auto_nomogram` applies, the axes take the
transforms it recovers. Every curve of `z` then comes out as a straight line, and the lines
are parallel. For a product, those transforms are logarithms, so you get log–log paper:

```@example conc
concurrency_chart((x, y) -> x * y, :x => (1, 10), :y => (1, 10), :z)
```

Straight lines are easier to draw and to read between, which is why charts drawn by hand
used to go to the trouble of finding such transforms.

## How the curves are chosen

Each field gets curves at round values across the range it takes on the chart, or across
its `domain` if it declares one. They are chosen the way an axis's ticks are, so `ticks`,
`tickstyle` and `tickformat` on a field's [`Variable`](@ref) all work.

A family's curves bunch together wherever the field changes quickly. Before a curve is
drawn at a point, it is checked against its neighbours: it is drawn only where it stands at
least [`NomogramStyle`](@ref)'s `min_curve_spacing` from the nearest curve that is at least
as round. The distance comes from the field's gradient, so it is known at every point
rather than estimated. A crowded family therefore thins out locally. The finer values stop
where they would run together, and the rounder ones carry on. The Moody diagram above shows
this on its left side.

The roundest curves are labelled where they leave the chart, past the right or top edge if
they can be. A label that would overlap another label, or cross a curve, is left off. Each
family's name goes beside its row of labels.

## Current limits

- Labels sit at the ends of curves. Charts that label curves along their length, inside
  the chart, aren't supported yet.
- Every family is drawn in the axis colour. Pass `gridcolor` to [`svg`](@ref) or
  [`nomogramplot`](@ref) to draw all the curves lighter; you can't yet colour families
  individually.
- The axes are straight and at right angles. Charts drawn on another grid, such as the
  Smith chart, aren't supported yet.
