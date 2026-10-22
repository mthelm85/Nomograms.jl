"""
    NomogramStyle(; kwargs...)

Sizes and thresholds for turning a [`Nomogram`](@ref) into geometry. All lengths are in
chart units -- the same units as a nomogram's `width` and `height`, which default to 70 by
100 -- so one style carries over unchanged between output sizes and backends.

# Keywords
- `majorticklength`, `minorticklength`: tick mark lengths.
- `ticklabelgap`: gap between a tick's outer end and its label.
- `ticklabelsize`, `titlesize`: text heights.
- `titlegap`: gap between the top of an axis and its title.
- `min_major_spacing`, `min_minor_spacing`: a tick is dropped if it would land closer
  than this to one already placed. This is what keeps a non-uniform axis legible.
- `label_minors`: whether minor ticks get labels too.
- `endpoint_ticks`: whether to add the ends of each domain as labelled ticks when they fit.
- `axislinewidth`, `tickwidth`, `isoplethwidth`: stroke widths.
- `isopleth_overhang`: how far an alignment line extends past the outermost scale.
- `markersize`: radius of the dot where an isopleth meets a scale.
"""
Base.@kwdef struct NomogramStyle
    majorticklength::Float64 = 3.0
    minorticklength::Float64 = 1.6
    ticklabelgap::Float64 = 1.2
    ticklabelsize::Float64 = 3.2
    titlesize::Float64 = 4.4
    titlegap::Float64 = 4.5
    min_major_spacing::Float64 = 4.0
    min_minor_spacing::Float64 = 1.3
    label_minors::Bool = false
    endpoint_ticks::Bool = true
    axislinewidth::Float64 = 0.6
    tickwidth::Float64 = 0.35
    isoplethwidth::Float64 = 0.5
    isopleth_overhang::Float64 = 6.0
    markersize::Float64 = 1.1
end

"""
    TickMark

One tick on a scale: where it sits, which way it points away from the axis, how long it is,
and what (if anything) it is labelled with. `level` is 1 for major ticks and 2 for minor.
"""
struct TickMark
    value::Float64
    position::Point
    normal::Point
    length::Float64
    label::Union{Nothing,String}
    level::Int
end

"""
    ScaleDrawing

Laid-out geometry for one scale: the axis polyline, its ticks, and where its title goes.
"""
struct ScaleDrawing
    name::Symbol
    label::String
    path::Vector{Point}
    ticks::Vector{TickMark}
    titlepos::Point
end

"""
    IsoplethDrawing

Laid-out geometry for one alignment line: the segment to stroke, and the points where it
meets each scale.
"""
struct IsoplethDrawing
    segment::Tuple{Point,Point}
    points::Vector{Point}
    label::String
end

"""
    NomogramDrawing

A complete, backend-independent description of a rendered nomogram: everything a renderer
needs, and nothing about how to draw it. [`svg`](@ref) and the Makie extension are both
thin walks over this structure, and tests can assert on tick positions without rendering
anything.
"""
struct NomogramDrawing
    scales::Vector{ScaleDrawing}
    isopleths::Vector{IsoplethDrawing}
    title::String
    titlepos::Point
    bbox::NTuple{4,Float64}
    style::NomogramStyle
end

_norm(p::Point) = hypot(p[1], p[2])
_dist(a::Point, b::Point) = _norm(_sub(a, b))

function _tangent(s::Scale, v::Float64)
    lo, hi = s.domain
    h = max(1e-6 * (hi - lo), 1e-8)
    a = max(v - h, lo)
    b = min(v + h, hi)
    b > a || return (0.0, 1.0)
    d = _sub(point(s, b), point(s, a))
    n = _norm(d)
    n == 0 && return (0.0, 1.0)
    return (d[1] / n, d[2] / n)
end

# Unit normal pointing to whichever side of the axis the ticks and labels belong on.
function _normal(s::Scale, v::Float64)
    t = _tangent(s, v)
    n1 = (t[2], -t[1])
    n2 = (-t[2], t[1])
    if s.side === :left
        return n1[1] <= n2[1] ? n1 : n2
    else
        return n1[1] >= n2[1] ? n1 : n2
    end
end

"""
    thin_ticks(s::Scale, majors, minors, style) -> Vector{TickMark}

Place as many of the candidate ticks as will fit, rejecting any that would land within
`min_major_spacing` / `min_minor_spacing` of one already placed.

Thinning happens in *page* distance, not value distance, which is what makes a reciprocal
or logarithmic axis readable: the tick values stay round numbers, and crowding near a
compressed end of the axis simply drops some of them.
"""
function thin_ticks(
    s::Scale,
    majors::AbstractVector{<:Real},
    minors::AbstractVector{<:Real},
    style::NomogramStyle,
)
    v = s.variable
    ticks = TickMark[]
    placed = Point[]

    function tryplace(value, level)
        value = float(value)
        p = point(s, value)
        gap = level == 1 ? style.min_major_spacing : style.min_minor_spacing
        all(_dist(p, q) >= gap for q in placed) || return false
        label = (level == 1 || style.label_minors) ? v.tickformat(value) : nothing
        len = level == 1 ? style.majorticklength : style.minorticklength
        push!(ticks, TickMark(value, p, _normal(s, value), len, label, level))
        push!(placed, p)
        return true
    end

    # The ends of the domain go down first: a reader needs to know where an axis stops,
    # so a round tick that crowds an endpoint is the one to sacrifice, not the endpoint.
    if style.endpoint_ticks
        for value in s.domain
            tryplace(value, 1)
        end
    end
    for value in majors
        tryplace(value, 1)
    end
    for value in minors
        tryplace(value, 2)
    end
    sort!(ticks; by = t -> t.value)
    return ticks
end

# Sample the scale curve, collapsing to a chord when it is straight -- which every axis of
# a parallel-layout nomogram is.
function _path(s::Scale; nsamples::Int = 129)
    lo, hi = s.domain
    a, b = point(s, lo), point(s, hi)
    pts = [point(s, x) for x in range(lo, hi; length = nsamples)]
    chord = _sub(b, a)
    len = _norm(chord)
    if len > 0
        tol = 1e-9 * max(len, 1.0)
        straight = all(abs(_cross(chord, _sub(p, a))) / len <= tol for p in pts)
        straight && return [a, b]
    end
    return pts
end

_textwidth(s::AbstractString, size::Real) = 0.6 * size * length(s)

function _expand!(bb::Vector{Float64}, p::Point; padx::Real = 0.0, pady::Real = 0.0)
    bb[1] = min(bb[1], p[1] - padx)
    bb[2] = min(bb[2], p[2] - pady)
    bb[3] = max(bb[3], p[1] + padx)
    bb[4] = max(bb[4], p[2] + pady)
    return bb
end

"""
    drawing(nom::Nomogram; style = NomogramStyle(), isopleths = Isopleth[])

Lay out a nomogram: choose and position ticks, build axis paths, and place titles and any
alignment lines. Returns a [`NomogramDrawing`](@ref) of pure geometry, shared by every
rendering backend.
"""
function drawing(
    nom::Nomogram;
    style::NomogramStyle = NomogramStyle(),
    isopleths::AbstractVector{Isopleth} = Isopleth[],
)
    bb = [Inf, Inf, -Inf, -Inf]
    scaledrawings = ScaleDrawing[]
    for s in nom.scales
        majors, minors = candidate_ticks(s.variable, s.domain)
        ticks = thin_ticks(s, majors, minors, style)
        path = _path(s)
        for p in path
            _expand!(bb, p)
        end
        for t in ticks
            outer = (
                t.position[1] + t.normal[1] * t.length,
                t.position[2] + t.normal[2] * t.length,
            )
            _expand!(bb, outer)
            isnothing(t.label) && continue
            reach = style.ticklabelgap + _textwidth(t.label, style.ticklabelsize)
            tip = (outer[1] + t.normal[1] * reach, outer[2] + t.normal[2] * reach)
            _expand!(bb, tip; pady = style.ticklabelsize / 2)
        end
        top = argmax(p -> p[2], path)
        titlepos = (top[1], top[2] + style.titlegap)
        _expand!(
            bb,
            (titlepos[1], titlepos[2] + style.titlesize);
            padx = _textwidth(s.variable.label, style.titlesize) / 2,
        )
        push!(scaledrawings, ScaleDrawing(name(s), s.variable.label, path, ticks, titlepos))
    end

    isodrawings = IsoplethDrawing[]
    for iso in isopleths
        pts = iso.points
        length(pts) >= 2 || continue
        a = argmin(p -> p[1], pts)
        b = argmax(p -> p[1], pts)
        d = _sub(b, a)
        n = _norm(d)
        u = n == 0 ? (1.0, 0.0) : (d[1] / n, d[2] / n)
        oh = style.isopleth_overhang
        seg = ((a[1] - u[1] * oh, a[2] - u[2] * oh), (b[1] + u[1] * oh, b[2] + u[2] * oh))
        _expand!(bb, seg[1])
        _expand!(bb, seg[2])
        if !isempty(iso.label)
            reach = style.ticklabelgap + _textwidth(iso.label, style.ticklabelsize)
            _expand!(bb, (seg[2][1] + reach, seg[2][2]); pady = style.ticklabelsize / 2)
        end
        push!(isodrawings, IsoplethDrawing(seg, pts, iso.label))
    end

    titlepos = (nom.width / 2, bb[4] + style.titlegap)
    if !isempty(nom.title)
        _expand!(
            bb,
            (titlepos[1], titlepos[2] + style.titlesize);
            padx = _textwidth(nom.title, style.titlesize) / 2,
        )
    end

    return NomogramDrawing(
        scaledrawings,
        isodrawings,
        nom.title,
        titlepos,
        (bb[1], bb[2], bb[3], bb[4]),
        style,
    )
end
