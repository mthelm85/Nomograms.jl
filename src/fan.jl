"""
    FanNomogram

A chart redrawn by [`fan`](@ref) with its axes as spokes from one centre. It wraps the
reshaped chart, which keeps its own type and does all the solving, and adds the rings drawn
around it.
"""
struct FanNomogram{N<:AbstractNomogram} <: AbstractNomogram
    nom::N
    scales::Vector{Scale}
    width::Float64
    height::Float64
    title::String
    rings::Vector{Vector{Point}}
end

_guides(f::FanNomogram) = f.rings

# Everything about reading the chart is the reshaped chart's own.
variablenames(f::FanNomogram) = variablenames(f.nom)
scale(f::FanNomogram, n::Symbol) = scale(f.nom, n)
evaluate(f::FanNomogram; values...) = evaluate(f.nom; values...)
evaluate_all(f::FanNomogram; values...) = evaluate_all(f.nom; values...)
isopleth(f::FanNomogram; kwargs...) = isopleth(f.nom; kwargs...)
isopleth(f::FanNomogram, solution::NamedTuple; kwargs...) =
    isopleth(f.nom, solution; kwargs...)
max_relation_residual(f::FanNomogram; kwargs...) = max_relation_residual(f.nom; kwargs...)

function Base.show(io::IO, f::FanNomogram)
    print(io, "FanNomogram(")
    show(io, f.nom)
    print(io, ")")
    return nothing
end

"""
    fan(nom; spread = 90, hub = 0.35, radius = 100, inward = :auto, rings = true) -> FanNomogram

Redraw a chart with its axes as spokes radiating from one centre, inside a pair of arcs --
the same chart, giving exactly the same answers, in a shape that suits a wall.

Parallel lines meet at a point at infinity, and a projective map can bring that point onto
the page. Collinear points stay collinear under such a map, which is why the chart still
works (see [`project`](@ref)); every vertical axis becomes a spoke through the new point,
the turning lines of a [`compound`](@ref) chart and the strands of a net included.

Straight lines stay straight, so the far ends of the spokes lie on a straight chord rather
than an arc, and the middle spokes stop short of the outer ring -- by more the wider the
fan. At 90° the ends line up level with the outer spokes' tips.

The map also squeezes whatever lies near the centre. So by default the chart is laid out
both ways round, with the top of each axis at the rim and at the hub, and whichever keeps
more tick labels is used: a scale that already crowds one end, such as a reciprocal one, is
best with that end at the rim.

# Keywords
- `spread = 90`: the angle in degrees between the outermost spokes; less than 180.
- `hub = 0.35`: where the spokes start, as a fraction of their length from the centre.
- `radius = 100`: the length of the outer spokes, in chart units.
- `inward = :auto`: `true` puts the top of each axis at the hub, `false` at the rim.
- `rings = true`: draw the arcs. They are [`NomogramDrawing`](@ref) guides, which labels
  keep clear of.
- `title`: defaults to the chart's own.

Works on any chart, but it is designed for parallel axes; a diagonal or curved scale is
reshaped too, and stays correct, but does not become a spoke.
"""
function fan(
    nom::AbstractNomogram;
    spread::Real = 90,
    hub::Real = 0.35,
    radius::Real = 100,
    inward::Union{Bool,Symbol} = :auto,
    rings::Bool = true,
    title::AbstractString = nom.title,
)
    0 < spread < 180 ||
        throw(ArgumentError("spread must be between 0 and 180 degrees, got $spread"))
    0 < hub < 1 || throw(ArgumentError("hub must be between 0 and 1, got $hub"))
    radius > 0 || throw(ArgumentError("radius must be positive, got $radius"))
    inward isa Bool ||
        inward === :auto ||
        throw(ArgumentError("inward must be true, false or :auto, got :$inward"))
    if inward === :auto
        out = fan(nom; spread, hub, radius, inward = false, rings, title)
        in_ = fan(nom; spread, hub, radius, inward = true, rings, title)
        labels(f) = count(
            t -> !isnothing(t.label),
            (t for sd in drawing(f).scales if sd.kind === :scale for t in sd.ticks),
        )
        return labels(in_) > labels(out) ? in_ : out
    end

    R, r = float(radius), float(hub) * radius
    h, θ0, θ1 = _fan_map(nom.scales, float(spread), r, R, inward)
    # Rings sit clear of the spokes' ends, where titles and end labels go -- far enough out
    # that a long title on a slanting outer spoke, which spreads sideways, still fits
    # under the ring rather than being pushed round to the hub.
    arc(ρ, a, b) = Point[(ρ * cosd(θ), ρ * sind(θ)) for θ in range(a, b; length = 97)]
    ringpaths =
        rings ?
        [arc(R + 24, θ0 + 4, θ1 - 4), arc(R + 25.5, θ0 + 4, θ1 - 4), arc(r - 9, θ0, θ1)] :
        Vector{Point}[]

    # Put the chart at the origin, as every construction does.
    pts = [h(point(s, v)) for s in nom.scales for v in range(s.domain...; length = 33)]
    append!(pts, (p for g in ringpaths for p in g))
    x0, y0 = minimum(first, pts), minimum(last, pts)
    x1, y1 = maximum(first, pts), maximum(last, pts)
    shift(p) = (p[1] - x0, p[2] - y0)
    g = shift ∘ h
    mapped = _remap(nom, g, x1 - x0, y1 - y0, String(title))
    return FanNomogram(
        mapped,
        mapped.scales,
        x1 - x0,
        y1 - y0,
        String(title),
        [shift.(path) for path in ringpaths],
    )
end

# The map: the left edge of the axes' extent to a ray at θ0, the right edge to one at θ1,
# and the bottom and top to radii r and R (or R and r, inward). Vertical lines then all
# pass through the origin.
function _fan_map(scales, spread, r, R, inward)
    pts = [point(s, v) for s in scales for v in range(s.domain...; length = 33)]
    x0, y0 = minimum(first, pts), minimum(last, pts)
    x1, y1 = maximum(first, pts), maximum(last, pts)
    (x1 > x0 && y1 > y0) ||
        throw(ArgumentError("a chart with no width or height cannot be fanned"))
    θ0, θ1 = 90 + spread / 2, 90 - spread / 2
    ray(ρ, θ) = (ρ * cosd(θ), ρ * sind(θ))
    lo, hi = inward ? (R, r) : (r, R)
    M = projective_transform(
        [(x0, y0), (x1, y0), (x1, y1), (x0, y1)],
        [ray(lo, θ0), ray(lo, θ1), ray(hi, θ1), ray(hi, θ0)],
    )
    h(p) = (w = M * [p[1], p[2], 1.0]; Point((w[1] / w[3], w[2] / w[3])))
    return h, θ0, θ1
end

_remap(s::Scale, g) = Scale(s.variable, s.domain, v -> g(s.curve(v)), s.side, s.kind)

_remap(nom::Nomogram, g, w, h, title) = Nomogram(
    [_remap(s, g) for s in nom.scales];
    width = w,
    height = h,
    title,
    relation = nom.relation,
)

function _remap(cn::CompoundNomogram, g, w, h, title)
    stages = [_remap(st, g, w, h, st.title) for st in cn.stages]
    return CompoundNomogram(stages, [_remap(s, g) for s in cn.scales], w, h, title)
end

function _remap(nn::NetNomogram, g, w, h, title)
    net = Net(nn.net.u, nn.net.v, (u, v) -> g(nn.net.surface(u, v)))
    return NetNomogram(
        net,
        _remap(nn.a, g),
        _remap(nn.b, g),
        [_remap(s, g) for s in nn.scales],
        w,
        h,
        title,
        nn.relation,
    )
end

function _remap(f::FanNomogram, g, w, h, title)
    throw(ArgumentError("this chart is already a fan; fan the original instead"))
end
