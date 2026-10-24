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
  variable of a construction, whose range is then derived from the other two. With
  Unitful loaded the ends may be quantities, `(1u"ft", 10u"ft")`, which also sets `unit`.
- `unit`: the unit the axis is drawn in, such as `u"ft"` (needs Unitful). Values passed
  to [`evaluate`](@ref) and [`isopleth`](@ref) may then be quantities in any unit of the
  same dimension, results come back as quantities, and a plain number is read as being in
  this unit. Defaults to the unit of the domain's ends, if they have one.
- `label`: axis title; defaults to `string(name)`. The unit is drawn after it,
  `length (ft)`, unless `unitlabel = false`.
- `ticks`, `minorticks`: explicit tick values, bypassing automatic selection.
- `tickformat`: `value -> String`. Defaults to [`fractionformat`](@ref) for
  `tickstyle = :fraction` and to plain numbers otherwise.
- `tickstyle`: how automatic tick values are picked.
  - `:linear` or `:log` for round decimals spaced evenly or by decade;
  - `:integer` for whole numbers only, for a count such as gear teeth or turns;
  - `:fraction` for fractions whose denominator divides `maxden` -- halves, quarters,
    eighths, sixteenths by default -- labelled as such: `1 3/8`, not `1.375`;
  - `:auto` to take it from the transform. A construction may override `:auto` when it
    knows better: an N-chart's ratio axis gets logarithmic ticks when it spans a decade or
    more.
  Whichever is chosen, ticks are thinned by distance on the page, roundest first.
- `maxden = 16`: for `tickstyle = :fraction`, the finest fraction offered; every
  denominator dividing it is used. 16 gives halves to sixteenths, 12 gives halves, thirds,
  quarters, sixths and twelfths.
- `side`: `:left`, `:right`, or `:auto` (the construction chooses).
- `titleend`: which end of the scale the title sits off — `:start` or `:stop` for the low
  or high end of the domain, or `:auto` for whichever is higher on the page. Worth setting
  where the automatic choice lands somewhere crowded, as it does for the strands of a
  [`Net`](@ref) that converge at one end.
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
    titleend::Symbol
    maxden::Int
    unit::Any
    unitlabel::Bool
end

# Units are handled by an extension, so the core only ever sees these four hooks. With
# no unit they do nothing: plain numbers in, plain numbers out.

# A value for a variable as a plain number in its display unit. A plain number is taken to
# be in that unit already; a quantity is converted, or refused if its dimension is wrong.
_strip(::Any, x::Real, name::Symbol) = float(x)
_strip(::Any, x, name::Symbol) =
    throw(ArgumentError("the value for :$name must be a number, got $(repr(x))"))
_attach(::Nothing, x) = x
_infer_unit(::Any) = nothing
_unitstring(::Nothing) = ""
_unitstring(u) = string(u)

_number(v::Variable, x) = _strip(v.unit, x, v.name)

"""
    axislabel(v::Variable) -> String

The title drawn for a variable's axis: its `label`, followed by its unit in brackets when
it has one and `unitlabel` is on.
"""
function axislabel(v::Variable)
    (isempty(v.label) || isnothing(v.unit) || !v.unitlabel) && return v.label
    return "$(v.label) ($(_unitstring(v.unit)))"
end

function Variable(
    name::Symbol,
    transform::ScaleTransform = LinearScale();
    domain = nothing,
    label::AbstractString = string(name),
    ticks = nothing,
    minorticks = nothing,
    tickformat::Union{Nothing,Function} = nothing,
    tickstyle::Symbol = :auto,
    side::Symbol = :auto,
    titleend::Symbol = :auto,
    maxden::Integer = 16,
    unit = nothing,
    unitlabel::Bool = true,
)
    side in (:left, :right, :auto) ||
        throw(ArgumentError("side must be :left, :right or :auto, got :$side"))
    tickstyle in (:auto, :linear, :log, :integer, :fraction) || throw(
        ArgumentError(
            "tickstyle must be :auto, :linear, :log, :integer or :fraction, got :$tickstyle",
        ),
    )
    maxden >= 1 || throw(ArgumentError("maxden must be at least 1, got $maxden"))
    if isnothing(tickformat)
        tickformat = tickstyle === :fraction ? fractionformat(maxden) : default_tickformat
    end
    titleend in (:auto, :start, :stop) ||
        throw(ArgumentError("titleend must be :auto, :start or :stop, got :$titleend"))
    if isnothing(unit) && !isnothing(domain)
        unit = something(_infer_unit(domain[1]), _infer_unit(domain[2]), Some(nothing))
    end
    num(x) = _strip(unit, x, name)
    dom = isnothing(domain) ? nothing : (num(domain[1]), num(domain[2]))
    if !isnothing(dom) && !(dom[1] < dom[2])
        throw(ArgumentError("domain of :$name must be increasing, got $dom"))
    end
    return Variable(
        name,
        transform,
        dom,
        String(label),
        isnothing(ticks) ? nothing : Float64[num(t) for t in ticks],
        isnothing(minorticks) ? nothing : Float64[num(t) for t in minorticks],
        tickformat,
        tickstyle,
        side,
        titleend,
        Int(maxden),
        unit,
        unitlabel,
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
    kind::Symbol
end

"""
    Scale(variable, domain, curve, side, kind = :scale)

`kind` is `:scale` for a scale a reader takes a value from, or `:grid` for one strand of a
[`Net`](@ref). Grid strands are drawn lighter, since a net is a backdrop rather than the
instrument.
"""
Scale(v::Variable, domain, curve, side::Symbol) = Scale(v, domain, curve, side, :scale)

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
            "]",
            isnothing(s.variable.unit) ? "" : " " * _unitstring(s.variable.unit),
            "  at x = ",
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
