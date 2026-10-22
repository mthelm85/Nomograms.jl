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

### Ticks

Tick values are chosen to be round numbers *in the variable*, then thinned by distance *on
the page* ([`thin_ticks`](@ref)). That split is what keeps a non-uniform axis readable: a
reciprocal or logarithmic scale keeps its labels where it is stretched and drops them
where it is compressed, and the numbers that survive are still round ones. Axis endpoints
and decade boundaries are placed first, so a crowded axis gives up `9` before it gives up
`10`.

Where an axis compresses so hard that automatic selection has nothing to work with, pass
`ticks = [...]` to the [`Variable`](@ref); anything that will not fit is still dropped.

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

## Reference

See the [Reference](@ref reference) page for the full API.
