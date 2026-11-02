"""
    fan(nom; spread = 90, hub = 0.35, radius = 100, inward = :auto) -> nomogram

Redraw a chart with its axes as spokes radiating from one centre. The fanned chart gives
exactly the same answers as the original. The result is the same kind of chart as `nom`: a
[`Nomogram`](@ref), a [`CompoundNomogram`](@ref) or a [`GridNomogram`](@ref).

Parallel lines meet at a point at infinity, and a projective map can bring that point onto
the page. Collinear points stay collinear under such a map, so the chart still works (see
[`project`](@ref)). Every vertical axis becomes a spoke through the new point, including the
pivot lines of a [`compound`](@ref) chart and the strands of a grid.

Straight lines stay straight, so the far ends of the spokes lie on a straight chord, not an
arc, and the middle spokes are shorter than the outer ones; the wider the fan, the bigger
the difference. At 90° the middle spokes' tips line up level with the outer spokes' tips.

The map also squeezes everything near the centre. So by default the chart is laid out both
ways round, once with the top of each axis at the rim and once with it at the hub, and the
layout that keeps more tick labels is used. A scale that is already crowded at one end, such
as a reciprocal scale, works best with that end at the rim.

# Keywords
- `spread = 90`: the angle in degrees between the outermost spokes; less than 180.
- `hub = 0.35`: where the spokes start, as a fraction of their length from the centre.
- `radius = 100`: the length of the outer spokes, in chart units.
- `inward = :auto`: `true` puts the top of each axis at the hub, `false` at the rim.
- `title`: defaults to the chart's own.

`fan` works on any chart, but it is designed for parallel axes. A diagonal or curved scale
is reshaped too and stays correct, but it does not become a spoke.
"""
function fan(
    nom::AbstractNomogram;
    spread::Real = 90,
    hub::Real = 0.35,
    radius::Real = 100,
    inward::Union{Bool,Symbol} = :auto,
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
        out = fan(nom; spread, hub, radius, inward = false, title)
        in_ = fan(nom; spread, hub, radius, inward = true, title)
        labels(f) = count(
            t -> !isnothing(t.label),
            (t for sd in drawing(f).scales if sd.kind === :scale for t in sd.ticks),
        )
        return labels(in_) > labels(out) ? in_ : out
    end

    R, r = float(radius), float(hub) * radius
    h = _fan_map(nom.scales, float(spread), r, R, inward)

    # Put the chart at the origin, as every construction does.
    pts = [h(point(s, v)) for s in nom.scales for v in range(s.domain...; length = 33)]
    x0, y0 = minimum(first, pts), minimum(last, pts)
    x1, y1 = maximum(first, pts), maximum(last, pts)
    g = p -> (q = h(p); Point((q[1] - x0, q[2] - y0)))
    return _remap(nom, g, x1 - x0, y1 - y0, String(title))
end

# The map: the left edge of the axes' extent to a ray at 90° + spread/2, the right edge to
# one at 90° - spread/2, and the bottom and top to radii r and R (or R and r, inward).
# Vertical lines then all pass through the origin.
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
    return p -> (w = M * [p[1], p[2], 1.0]; Point((w[1] / w[3], w[2] / w[3])))
end

fan(c::ContourChart; kwargs...) = throw(
    ArgumentError(
        "`fan` turns the parallel axes of an alignment chart into spokes. A contour chart " *
        "has no such axes; to draw one on a disc, use `axes = false` and `rings`.",
    ),
)

# `project` on a chart of several parts: the same map, and the same fit to the canvas, for
# every stage, strand and surface, so that the parts stay aligned with one another. A
# homography keeps straight lines straight, so each stage still computes what it did.
function project(
    nom::Union{CompoundNomogram,GridNomogram},
    H::AbstractMatrix;
    width::Real = nom.width,
    height::Real = nom.height,
    title::AbstractString = nom.title,
    preserve_aspect::Bool = false,
    nsamples::Int = 129,
)
    _check_homography(H, width, height)
    W, Hc = float(width), float(height)
    pts = [
        _project_point(H, point(s, v), name(s)) for s in nom.scales for
        v in range(s.domain...; length = nsamples)
    ]
    x0, x1 = extrema(first, pts)
    y0, y1 = extrema(last, pts)
    (x1 > x0 && y1 > y0) ||
        throw(ArgumentError("the projective map collapses the chart onto a line"))
    ax, ay = W / (x1 - x0), Hc / (y1 - y0)
    preserve_aspect && (ax = ay = min(ax, ay))
    bx = (W - ax * (x1 - x0)) / 2 - ax * x0
    by = (Hc - ay * (y1 - y0)) / 2 - ay * y0
    g = function (p)
        q = _project_point(H, (float(p[1]), float(p[2])), :chart)
        return Point((ax * q[1] + bx, ay * q[2] + by))
    end
    return _remap(nom, g, W, Hc, String(title))
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

function _remap(nn::GridNomogram, g, w, h, title)
    grid = Grid(nn.grid.u, nn.grid.v, (u, v) -> g(nn.grid.surface(u, v)))
    return GridNomogram(
        grid,
        _remap(nn.a, g),
        _remap(nn.b, g),
        [_remap(s, g) for s in nn.scales],
        w,
        h,
        title,
        nn.relation,
    )
end
