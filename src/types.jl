"""
    Point

A position in chart coordinates, `(x, y)` with `y` increasing upwards.
"""
const Point = NTuple{2,Float64}

"""
    Variable(name, transform = LinearScale(); kwargs...)

One variable of a nomogram, before it has been placed on the page.

# Arguments
- `name::Symbol`: the keyword used to refer to this variable in [`evaluate`](@ref) and
  [`isopleth`](@ref).
- `transform::ScaleTransform`: the `f` that this variable contributes to the separable
  relation `f1(u) + f2(v) = f3(w)`.

# Keywords
- `domain`: `(lo, hi)` in the variable's own units. May be `nothing` for the *result*
  variable of a construction, whose range is then derived from the other two.
- `label`: axis title; defaults to `string(name)`.
- `ticks`, `minorticks`: explicit tick values, bypassing automatic selection.
- `tickformat`: `value -> String`.
- `tickstyle`: `:linear` or `:log` to choose how automatic tick values are picked, or
  `:auto` to take it from the transform. A construction may override `:auto` when it knows
  better: an N-chart's ratio axis gets logarithmic ticks when it spans a decade or more.
- `side`: `:left`, `:right`, or `:auto` (the construction chooses).
"""
struct Variable
    name::Symbol
    transform::ScaleTransform
    domain::Union{Nothing,Tuple{Float64,Float64}}
    label::String
    ticks::Union{Nothing,Vector{Float64}}
    minorticks::Union{Nothing,Vector{Float64}}
    tickformat::Function
    tickstyle::Symbol
    side::Symbol
end

function Variable(
    name::Symbol,
    transform::ScaleTransform = LinearScale();
    domain = nothing,
    label::AbstractString = string(name),
    ticks = nothing,
    minorticks = nothing,
    tickformat::Function = default_tickformat,
    tickstyle::Symbol = :auto,
    side::Symbol = :auto,
)
    side in (:left, :right, :auto) ||
        throw(ArgumentError("side must be :left, :right or :auto, got :$side"))
    tickstyle in (:auto, :linear, :log) ||
        throw(ArgumentError("tickstyle must be :auto, :linear or :log, got :$tickstyle"))
    dom = isnothing(domain) ? nothing : (float(domain[1]), float(domain[2]))
    if !isnothing(dom) && !(dom[1] < dom[2])
        throw(ArgumentError("domain of :$name must be increasing, got $dom"))
    end
    return Variable(
        name,
        transform,
        dom,
        String(label),
        isnothing(ticks) ? nothing : collect(float.(ticks)),
        isnothing(minorticks) ? nothing : collect(float.(minorticks)),
        tickformat,
        tickstyle,
        side,
    )
end

"""
    resolved_tickstyle(v::Variable) -> Symbol

`:linear` or `:log`, taking the variable's own `tickstyle` when it is not `:auto` and
otherwise asking the transform.
"""
resolved_tickstyle(v::Variable) =
    v.tickstyle === :auto ? tickstyle(v.transform) : v.tickstyle

"""
    Scale

A [`Variable`](@ref) placed on the page: a parametric curve `value -> (x, y)` over a
resolved domain. Collinearity of one point per scale is the whole content of a nomogram, so
every layout — parallel, curved, or hand-built — reduces to a set of these.
"""
struct Scale
    variable::Variable
    domain::Tuple{Float64,Float64}
    curve::Function
    side::Symbol
end

"""
    point(s::Scale, value) -> Point

Position of `value` on scale `s`, in chart coordinates.
"""
point(s::Scale, value) = Point(s.curve(float(value)))

name(s::Scale) = s.variable.name
domain(s::Scale) = s.domain

"""
    AbstractNomogram

A chart made of [`Scale`](@ref)s on a canvas. Everything about *drawing* one -- ticks,
labels, bounding box, both backends -- needs only the scales, the canvas and the title, so
it is written against this type. What differs between subtypes is how they are *solved*: a
[`Nomogram`](@ref) is one alignment line through every scale, while a
[`CompoundNomogram`](@ref) is several stages read in turn.
"""
abstract type AbstractNomogram end

"""
    Nomogram

A set of [`Scale`](@ref)s on a shared canvas, plus the relation they encode. Every scale is
crossed by one alignment line.

Build one with [`additive_nomogram`](@ref) or [`multiplicative_nomogram`](@ref), or
directly from scales for layouts those constructions do not cover. For a relation needing
more than three variables, see [`compound`](@ref).
"""
struct Nomogram <: AbstractNomogram
    scales::Vector{Scale}
    width::Float64
    height::Float64
    title::String
    relation::Union{Nothing,Function}
end

function Nomogram(
    scales::AbstractVector{Scale};
    width::Real,
    height::Real,
    title::AbstractString = "",
    relation = nothing,
)
    length(scales) >= 3 ||
        throw(ArgumentError("a nomogram needs at least 3 scales, got $(length(scales))"))
    names = [name(s) for s in scales]
    allunique(names) || throw(ArgumentError("scale names must be unique, got $names"))
    return Nomogram(collect(scales), float(width), float(height), String(title), relation)
end

Base.getindex(nom::AbstractNomogram, n::Symbol) = scale(nom, n)

"""
    scale(nom::AbstractNomogram, name::Symbol) -> Scale

The scale carrying the given variable name.
"""
function scale(nom::AbstractNomogram, n::Symbol)
    i = findfirst(s -> name(s) === n, nom.scales)
    isnothing(i) && throw(
        KeyError(
            "no variable :$n in this nomogram (have " *
            join(":" .* string.(variablenames(nom)), ", ") *
            ")",
        ),
    )
    return nom.scales[i]
end

"""
    variablenames(nom::AbstractNomogram) -> Vector{Symbol}
"""
variablenames(nom::AbstractNomogram) = [name(s) for s in nom.scales]

"""
    Isopleth

A straight line laid across a [`Nomogram`](@ref) through one point per scale: the visual
form of a single solution of the relation. Build one with [`isopleth`](@ref).

# Fields
- `values::Dict{Symbol,Float64}`: the value read off each scale.
- `points::Vector{Point}`: where the line meets each scale, in scale order.
- `label::String`
"""
struct Isopleth
    values::Dict{Symbol,Float64}
    points::Vector{Point}
    label::String
    group::Int
end

Isopleth(values, points, label) = Isopleth(values, points, String(label), 0)

Base.getindex(iso::Isopleth, n::Symbol) = iso.values[n]

function Base.show(io::IO, nom::Nomogram)
    print(io, "Nomogram(")
    join(io, (":$(name(s))" for s in nom.scales), ", ")
    print(io, ")")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", nom::AbstractNomogram)
    println(
        io,
        "Nomogram: $(length(nom.scales)) scales, ",
        "canvas $(nom.width) x $(nom.height)",
    )
    isempty(nom.title) || println(io, "  title: ", nom.title)
    for s in nom.scales
        lo, hi = s.domain
        p = point(s, lo)
        println(
            io,
            "  :",
            rpad(string(name(s)), 10),
            " domain [",
            default_tickformat(lo),
            ", ",
            default_tickformat(hi),
            "]  at x = ",
            default_tickformat(p[1]),
        )
    end
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", iso::Isopleth)
    println(io, "Isopleth", isempty(iso.label) ? "" : " \"$(iso.label)\"", ":")
    for (k, v) in sort(collect(iso.values); by = first)
        println(io, "  :", rpad(string(k), 10), " = ", default_tickformat(v))
    end
    return nothing
end
