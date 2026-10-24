"""
    Nomograms

Build, solve and render nomograms (alignment charts).

A nomogram is an analog computer printed on paper. Lay a straightedge through the values
you know, and it crosses the remaining scale at the answer — no power, no keystrokes, just
geometry. It works because its three or more scales are drawn so that any straight line
crossing them picks out a solution of an equation. The entire geometric content is
d'Ocagne's collinearity criterion,

```
| x1  y1  1 |
| x2  y2  1 | = 0
| x3  y3  1 |
```

so this package represents every scale as a parametric curve `value -> (x, y)` and does all
its arithmetic by asking where a line crosses a curve. One solver then serves forward and
inverse problems alike, and layouts beyond the built-in parallel-axis construction need no
new machinery.

# Getting started

```julia
using Nomograms

nom = additive_nomogram(
    Variable(:R1, ReciprocalScale(); domain = (1, 10), label = "R1 (ohm)"),
    Variable(:R2, ReciprocalScale(); domain = (1, 10), label = "R2 (ohm)"),
    Variable(:R, ReciprocalScale(); label = "R (ohm)");
    title = "Resistors in parallel",
)

evaluate(nom; R1 = 4, R2 = 6)          # -> (R1 = 4.0, R2 = 6.0, R = 2.4)
save_svg("parallel.svg", nom; isopleths = [isopleth(nom; R1 = 4, R2 = 6)])
```

See [`additive_nomogram`](@ref), [`multiplicative_nomogram`](@ref), [`evaluate`](@ref),
[`isopleth`](@ref), [`svg`](@ref) and [`nomogramplot`](@ref).
"""
module Nomograms

using Printf: @sprintf

using FunctionWrappers: FunctionWrapper
import LinearAlgebra
import PrecompileTools
import Roots

export
    # transforms
    ScaleTransform,
    LinearScale,
    LogScale,
    ReciprocalScale,
    PowerScale,
    FunctionScale,
    ComposedScale,
    transform,
    invtransform,
    # model
    Variable,
    Scale,
    AbstractNomogram,
    Nomogram,
    CompoundNomogram,
    Net,
    NetNomogram,
    Isopleth,
    point,
    variablenames,
    # construction
    additive_nomogram,
    multiplicative_nomogram,
    ratio_nomogram,
    determinant_nomogram,
    project,
    projective_transform,
    fan,
    FanNomogram,
    compound,
    net_nomogram,
    netpoint,
    auto_nomogram,
    separability_defect,
    # solving
    evaluate,
    evaluate_all,
    isopleth,
    align_residual,
    solve_on,
    solve_all_on,
    relation_residual,
    max_relation_residual,
    # rendering
    NomogramStyle,
    NomogramDrawing,
    drawing,
    svg,
    save_svg,
    nomogramplot,
    nomogramplot!

include("transforms.jl")
include("types.jl")
include("ticks.jl")
include("construct.jl")
include("solve.jl")
include("compound.jl")
include("net.jl")
include("auto.jl")
include("auto_chain.jl")
include("drawing.jl")
include("svg.jl")
include("fan.jl")
include("precompile.jl")

"""
    nomogramplot(nom::Nomogram; kwargs...) -> Figure
    nomogramplot!(ax, nom::Nomogram; kwargs...) -> ax

Draw a nomogram with Makie. Available once a Makie backend is loaded:

```julia
using Nomograms, CairoMakie
fig = nomogramplot(nom; isopleths = [isopleth(nom; R1 = 4, R2 = 6)])
save("chart.pdf", fig)
```

# Keywords
- `isopleths`, `style`: as for [`drawing`](@ref).
- `axiscolor`, `textcolor`, `isoplethcolors`: colours.
- `pxperunit = 4.0`: figure size per chart unit (`nomogramplot` only).

For vector output with no plotting dependency at all, use [`svg`](@ref) instead.
"""
function nomogramplot(args...; kwargs...)
    throw(
        ArgumentError(
            "nomogramplot requires Makie: run `using CairoMakie` (or another " *
            "Makie backend) alongside Nomograms to load the plotting extension. " *
            "For dependency-free vector output, use `svg` or `save_svg`.",
        ),
    )
end

"""
    nomogramplot!(ax, nom::Nomogram; kwargs...)

Draw a nomogram into an existing Makie axis. See [`nomogramplot`](@ref).
"""
function nomogramplot!(args...; kwargs...)
    return nomogramplot(args...; kwargs...)
end

end # module
