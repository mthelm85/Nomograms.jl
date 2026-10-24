"""
    NomogramStyle(; kwargs...)

Sizes and thresholds for turning a [`Nomogram`](@ref) into geometry. All lengths are in
chart units -- the same units as a nomogram's `width` and `height`, which default to 70 by
100 -- so one style carries over unchanged between output sizes and backends.

The defaults follow the proportions of a classic d3 axis: hairline rules, short ticks and
small type. At the default `pxperunit = 4` a chart unit is 4 pixels, so they come out as
roughly 1px strokes, 6px major ticks, 3px of label padding and 10px labels.

# Keywords
- `majorticklength`, `minorticklength`: tick mark lengths.
- `graded_ticks`: whether ticks between those lengths are graded by importance, as on a
  ruler -- wholes longest, then halves, quarters and so on on a fraction axis; 10s, 5s,
  then the rest on an integer one; decades, whole digits, tenths on a logarithmic one.
  Labels sit the same distance from the axis whatever their tick's length. `false` gives
  every labelled tick the major length and every other the minor.
- `ticklabelgap`: gap between a tick's outer end and its label.
- `ticklabelsize`, `titlesize`: text heights.
- `titlegap`: gap between the end of an axis and its title.
- `min_major_spacing`, `min_minor_spacing`: a tick is dropped if it would land closer
  than this to one already placed. This is what keeps a non-uniform axis legible.
- `label_minors`: whether minor ticks get labels too.
- `endpoint_ticks`: whether to add the ends of each domain as labelled ticks when they fit.
  An end needing more than three significant figures is left to the axis line: that is a
  derived reach such as `237.4`, not a value anyone chose.
- `ticklabelrotation`: `:horizontal` leaves tick labels upright; `:tangent` turns each to
  run along its scale, which reads better on a steep or curved axis.
- `resolve_collisions`: whether to drop tick labels that would overlap a label from
  another scale, and nudge isopleth labels clear of both.
- `avoid_isopleths`: whether a tick label struck through by an alignment line is dropped
  too. A label a line runs across cannot be read, and the tick mark still shows where the
  value falls.
- `collisionpad`: extra clearance required between two labels.
- `axislinewidth`, `gridlinewidth`, `tickwidth`, `isoplethwidth`: stroke widths. The
  grid width applies to the strands of a [`Net`](@ref), which are drawn lighter than the
  scales a reader takes values from.
- `isopleth_overhang`: how far an alignment line extends past the outermost scale.
- `markersize`: radius of the dot where an isopleth meets a scale.
"""
Base.@kwdef struct NomogramStyle
    majorticklength::Float64 = 1.5
    minorticklength::Float64 = 0.8
    graded_ticks::Bool = true
    ticklabelgap::Float64 = 0.8
    ticklabelsize::Float64 = 2.5
    titlesize::Float64 = 2.9
    titlegap::Float64 = 3.0
    min_major_spacing::Float64 = 3.2
    min_minor_spacing::Float64 = 1.0
    label_minors::Bool = false
    endpoint_ticks::Bool = true
    ticklabelrotation::Symbol = :horizontal
    resolve_collisions::Bool = true
    avoid_isopleths::Bool = true
    collisionpad::Float64 = 0.3
    axislinewidth::Float64 = 0.25
    gridlinewidth::Float64 = 0.15
    tickwidth::Float64 = 0.25
    isoplethwidth::Float64 = 0.4
    isopleth_overhang::Float64 = 5.0
    markersize::Float64 = 0.8
end

"""
    TickMark

One tick on a scale: where it sits, which way it points away from the axis, how long it is,
and what (if anything) it is labelled with. `level` is 1 for major ticks and 2 for minor.

`labelpos`, `align` and `rotation` fully determine where the text goes, so every backend
places it identically and the layout can reason about the box it occupies.
"""
struct TickMark
    value::Float64
    position::Point
    normal::Point
    length::Float64
    label::Union{Nothing,String}
    level::Int
    labelpos::Point
    align::NTuple{2,Symbol}
    rotation::Float64
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
    titlealign::NTuple{2,Symbol}
    kind::Symbol
end

"""
    IsoplethDrawing

Laid-out geometry for one alignment line: the segment to stroke, the points where it meets
each scale, and where its label sits.
"""
struct IsoplethDrawing
    segment::Tuple{Point,Point}
    points::Vector{Point}
    label::String
    labelpos::Point
    labelalign::NTuple{2,Symbol}
    group::Int
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
    titlealign::NTuple{2,Symbol}
    bbox::NTuple{4,Float64}
    style::NomogramStyle
end

_segments(path::Vector{Point}) = [(path[k], path[k+1]) for k = 1:(length(path)-1)]

_norm(p::Point) = hypot(p[1], p[2])
_dist(a::Point, b::Point) = _norm(_sub(a, b))
_unit(p::Point) = (n = _norm(p); n == 0 ? (1.0, 0.0) : (p[1] / n, p[2] / n))
_offset(p::Point, d::Point, k::Real) = (p[1] + d[1] * k, p[2] + d[2] * k)

function _tangent(s::Scale, v::Float64)
    lo, hi = s.domain
    h = max(1e-6 * (hi - lo), 1e-8)
    a = max(v - h, lo)
    b = min(v + h, hi)
    b > a || return (0.0, 1.0)
    return _unit(_sub(point(s, b), point(s, a)))
end

# The two candidate normals are negatives of each other, so ordering them lexicographically
# by (x, y) names them unambiguously: `:left` is the smaller, `:right` the larger. Ordering
# on x alone leaves a horizontal scale with no way to say "label above" -- both sides would
# resolve to the same normal.
_normal_lt(a::Point, b::Point) = a[1] < b[1] || (a[1] == b[1] && a[2] < b[2])

"""
    side_normals(tangent) -> (left, right)

The two unit normals to a scale, named as [`Variable`](@ref)'s `side` names them.
"""
function side_normals(t::Point)
    n1 = (t[2], -t[1])
    n2 = (-t[2], t[1])
    return _normal_lt(n1, n2) ? (n1, n2) : (n2, n1)
end

# Unit normal pointing to whichever side of the axis the ticks and labels belong on.
function _normal(s::Scale, v::Float64)
    left, right = side_normals(_tangent(s, v))
    return s.side === :left ? left : right
end

"""
    align_for(direction) -> (halign, valign)

Which corner of a text box sits on its anchor point, given the direction the text is
pushed away from whatever it labels.

Both components matter. Using only the horizontal one -- enough while every axis is
vertical -- centres a label vertically on its anchor, so a label pushed straight down from
a tick on a curved scale would sit half on top of it.
"""
function align_for(d::Point)
    h = d[1] < -0.3 ? :right : d[1] > 0.3 ? :left : :center
    v = d[2] < -0.3 ? :top : d[2] > 0.3 ? :bottom : :center
    return (h, v)
end

# Text set along the scale rather than upright. The baseline is chosen so that the
# rotated frame's "up" is the outward normal, then flipped if that would read backwards.
function _tangent_label(n::Point)
    t = (n[2], -n[1])
    valign = :bottom
    if t[1] < -1e-12 || (abs(t[1]) <= 1e-12 && t[2] < 0)
        t = (-t[1], -t[2])
        valign = :top
    end
    return (atan(t[2], t[1]), (:center, valign))
end

_textwidth(s::AbstractString, size::Real) = 0.6 * size * length(s)

"""
    label_box(pos, text, size, align, rotation) -> (xmin, ymin, xmax, ymax)

The axis-aligned box a piece of text occupies, used both to size the canvas and to detect
labels running into each other.
"""
function label_box(
    pos::Point,
    text::AbstractString,
    size::Real,
    align::NTuple{2,Symbol},
    rotation::Real,
)
    w = _textwidth(text, size)
    h = float(size)
    x0 = align[1] === :left ? 0.0 : align[1] === :right ? -w : -w / 2
    y0 = align[2] === :bottom ? 0.0 : align[2] === :top ? -h : -h / 2
    corners = ((x0, y0), (x0 + w, y0), (x0, y0 + h), (x0 + w, y0 + h))
    c, s = cos(rotation), sin(rotation)
    xs = ntuple(i -> pos[1] + c * corners[i][1] - s * corners[i][2], 4)
    ys = ntuple(i -> pos[2] + s * corners[i][1] + c * corners[i][2], 4)
    return (minimum(xs), minimum(ys), maximum(xs), maximum(ys))
end

"""
    segment_hits_box(a, b, box, pad) -> Bool

Whether the segment from `a` to `b` passes through `box`, grown by `pad`.

Liang--Barsky: clip the segment against each of the box's four slabs in turn, and see
whether any of it is left. Label collisions are not only label against label -- an
alignment line drawn straight through a number makes it just as unreadable.
"""
function segment_hits_box(a::Point, b::Point, box::NTuple{4,Float64}, pad::Real)
    xmin, ymin = box[1] - pad, box[2] - pad
    xmax, ymax = box[3] + pad, box[4] + pad
    t0, t1 = 0.0, 1.0
    dx, dy = b[1] - a[1], b[2] - a[2]
    for (p, q) in
        ((-dx, a[1] - xmin), (dx, xmax - a[1]), (-dy, a[2] - ymin), (dy, ymax - a[2]))
        if p == 0
            q < 0 && return false      # parallel to this slab and outside it
        else
            r = q / p
            if p < 0
                r > t1 && return false
                r > t0 && (t0 = r)
            else
                r < t0 && return false
                r < t1 && (t1 = r)
            end
        end
    end
    return t0 <= t1
end

_overlaps(a::NTuple{4,Float64}, b::NTuple{4,Float64}, pad::Real) =
    a[1] - pad < b[3] && b[1] - pad < a[3] && a[2] - pad < b[4] && b[2] - pad < a[4]

"""
    thin_ticks(s::Scale, majors, minors, style) -> Vector{TickMark}

Place as many of the candidate ticks as will fit, rejecting any that would land within
`min_major_spacing` / `min_minor_spacing` of one already placed.

Thinning happens in *page* distance, not value distance, which is what makes a reciprocal
or logarithmic axis readable: the tick values stay round numbers, and crowding near a
compressed end of the axis simply drops some of them.

This only ever compares a scale against itself; labels colliding with *another* scale are
handled later, over the whole chart at once.
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

    function tryplace(value, level; endpoint = false)
        value = float(value)
        p = point(s, value)
        gap = level == 1 ? style.min_major_spacing : style.min_minor_spacing
        all(_dist(p, q) >= gap for q in placed) || return false
        label = (level == 1 || style.label_minors) ? v.tickformat(value) : nothing
        maj, mnr = style.majorticklength, style.minorticklength
        len = if !style.graded_ticks || endpoint
            # An axis end marks where the axis stops, so it is drawn in full.
            (level == 1 || endpoint) ? maj : mnr
        else
            tier, deepest = _tick_tier(v, s.domain, value, level)
            deepest == 0 ? maj : maj - (maj - mnr) * tier / deepest
        end
        n = _normal(s, value)
        # Labels keep one distance from the axis, so a column of them stays straight
        # however their ticks are graded.
        labelpos =
            _offset(p, n, (isnothing(label) ? len : max(len, maj)) + style.ticklabelgap)
        rotation, align = if style.ticklabelrotation === :tangent
            _tangent_label(n)
        else
            (0.0, align_for(n))
        end
        push!(ticks, TickMark(value, p, n, len, label, level, labelpos, align, rotation))
        push!(placed, p)
        return true
    end

    # The ends of the domain go down first: a reader needs to know where an axis stops,
    # so a round tick that crowds an endpoint is the one to sacrifice, not the endpoint.
    # An explicitly empty tick list is a request for a bare line, though -- a compound
    # nomogram's turning line is usually drawn blank -- so it overrides this.
    blank = !isnothing(v.ticks) && isempty(v.ticks)
    if style.endpoint_ticks && !blank
        for value in s.domain
            _readable_end(v, value) && tryplace(value, 1; endpoint = true)
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

# Whether an end of a domain is worth a label. A declared end is a round number someone
# chose; a derived one -- the reach of a result axis -- is wherever the arithmetic landed,
# and `237.4` or `1.883e6` at the end of an axis is noise rather than a reading. Three
# significant figures keeps every end anyone is likely to type, and `1.01` with it.
function _readable_end(x::Real)
    (iszero(x) || !isfinite(x)) && return true
    return isapprox(round(x; sigdigits = 3), x; rtol = 1e-9)
end

# On an integer or fraction axis the test is the style's own: 1 5/16 is a fine end for a
# ruler in sixteenths though it takes five figures as a decimal, and 2.5 is no end for a
# count however round it looks.
function _readable_end(v::Variable, x::Real)
    v.tickstyle in (:integer, :fraction) && return _style_allows(v, x)
    return _readable_end(x)
end

_unlabel(t::TickMark) = TickMark(
    t.value,
    t.position,
    t.normal,
    t.length,
    nothing,
    t.level,
    t.labelpos,
    t.align,
    t.rotation,
)

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

# Where a scale's title can go: off an end of the axis, continuing in the direction the
# curve was already heading. For a vertical axis that is straight up from the top, which
# is what it has always been; for a diagonal or a curve it follows the scale instead of
# jumping to whatever point happens to be highest. Both ends are offered, preferred first,
# unless `titleend` names one.
function _title_ends(s::Scale, path::Vector{Point}, style::NomogramStyle)
    ends = (
        (path[1], _sub(path[1], path[min(2, end)])),
        (path[end], _sub(path[end], path[max(end - 1, 1)])),
    )
    # Prefer the higher end, breaking ties on the left so the choice is deterministic.
    first =
        ends[1][1][2] > ends[2][1][2] + 1e-12 ? 1 :
        ends[2][1][2] > ends[1][1][2] + 1e-12 ? 2 : ends[1][1][1] <= ends[2][1][1] ? 1 : 2
    order =
        s.variable.titleend === :start ? (1,) :
        s.variable.titleend === :stop ? (2,) : (first, 3 - first)
    return map(order) do k
        at, dir = ends[k]
        away = _unit(dir)
        (_offset(at, away, style.titlegap), align_for(away), away)
    end
end

# Where a scale's title goes when nothing is in the way: its preferred end.
function _title_anchor(s::Scale, path::Vector{Point}, style::NomogramStyle)
    pos, align, away = first(_title_ends(s, path, style))
    return (pos, align, away)
end

function _expand_box!(bb::Vector{Float64}, box::NTuple{4,Float64})
    bb[1] = min(bb[1], box[1])
    bb[2] = min(bb[2], box[2])
    bb[3] = max(bb[3], box[3])
    bb[4] = max(bb[4], box[4])
    return bb
end

function _expand!(bb::Vector{Float64}, p::Point)
    return _expand_box!(bb, (p[1], p[2], p[1], p[2]))
end

"""
    nudge_clear(text, anchor, dir, size, align, style, occupied) -> (position, box)

Place a label at `anchor`, or at the nearest offset along `±dir` that clears everything in
`occupied`, and record the box it takes.

This is for text that must not simply be dropped: an axis title, or an isopleth label the
caller asked for by name. Tick labels are the other case -- there are plenty of those, so
one that will not fit loses its label and keeps its tick.
"""
function nudge_clear(
    text::AbstractString,
    anchor::Point,
    dir::Point,
    size::Real,
    align::NTuple{2,Symbol},
    style::NomogramStyle,
    occupied::Vector{NTuple{4,Float64}},
)
    if !style.resolve_collisions
        box = label_box(anchor, text, size, align, 0.0)
        push!(occupied, box)
        return (anchor, box)
    end
    step = float(size) * 1.25
    # Far enough out along `dir` is clear of everything, so this always ends with a place
    # that overlaps nothing -- never with the least bad overlap.
    kmax = _steps_clear(label_box(anchor, text, size, align, 0.0), occupied, step)
    for k = 0:kmax, sgn in (k == 0 ? (1,) : (1, -1))
        pos = _offset(anchor, dir, sgn * step * k)
        box = label_box(pos, text, size, align, 0.0)
        if !any(b -> _overlaps(box, b, style.collisionpad), occupied)
            push!(occupied, box)
            return (pos, box)
        end
    end
    error("unreachable: a label $(kmax) steps out still overlapped something")
end

# How many steps of `step` move `box` clear of everything in `boxes` and `points`,
# whichever way it goes: enough to cross the whole extent of both.
function _steps_clear(box, boxes, step; points = Point[])
    xs = Float64[box[1], box[3]]
    ys = Float64[box[2], box[4]]
    for b in boxes
        push!(xs, b[1], b[3])
        push!(ys, b[2], b[4])
    end
    for p in points
        push!(xs, p[1])
        push!(ys, p[2])
    end
    extent = hypot(maximum(xs) - minimum(xs), maximum(ys) - minimum(ys))
    return ceil(Int, extent / step) + 2
end

# Titles go down before anything else, against each other and against the axis lines of
# every other scale a reader takes values from. Each tries its preferred end, then the
# other end, then rows further out from either, until one is clear -- which far enough out
# is guaranteed, so a title is never drawn overlapping another or across an axis. Moving
# outward along the axis keeps a vertical axis's title centred over it, where sliding
# sideways would leave a wide title closer to a neighbour's axis than its own.
function _place_titles(
    nom::AbstractNomogram,
    paths::Vector{Vector{Point}},
    style::NomogramStyle,
    occupied::Vector{NTuple{4,Float64}},
)
    walls = [
        (i, seg) for
        (i, (s, path)) in enumerate(zip(nom.scales, paths)) if s.kind === :scale for
        seg in _segments(path)
    ]
    placed = Tuple{Point,NTuple{2,Symbol}}[]
    for (i, s) in enumerate(nom.scales)
        ends = _title_ends(s, paths[i], style)
        text = axislabel(s.variable)
        if isempty(text) || !style.resolve_collisions
            pos, align, _ = ends[1]
            isempty(text) ||
                push!(occupied, label_box(pos, text, style.titlesize, align, 0.0))
            push!(placed, (pos, align))
            continue
        end
        clear(box) =
            !any(b -> _overlaps(box, b, style.collisionpad), occupied) &&
            !any(
                w ->
                    w[1] != i &&
                    segment_hits_box(w[2][1], w[2][2], box, style.collisionpad),
                walls,
            )
        step = style.titlesize * 1.3
        kmax = _steps_clear(
            label_box(ends[1][1], text, style.titlesize, ends[1][2], 0.0),
            occupied,
            step;
            points = [p for (_, seg) in walls for p in seg],
        )
        spot = nothing
        for k = 0:kmax, (anchor, align, away) in ends
            pos = _offset(anchor, away, k * step)
            box = label_box(pos, text, style.titlesize, align, 0.0)
            if clear(box)
                spot = (pos, align, box)
                break
            end
        end
        isnothing(spot) && error("unreachable: title :$(name(s)) found no clear place")
        push!(occupied, spot[3])
        push!(placed, (spot[1], spot[2]))
    end
    return placed
end

# Tick labels are thinned per scale by distance along their own axis, which says nothing
# about labels belonging to other scales. Non-parallel layouts bring scales together, so
# compare the actual text boxes across the whole chart: titles come first, then major
# labels, then minor ones, and a label that cannot fit is dropped while its tick mark
# stays. `occupied` already holds the titles, which were placed clear of each other. A
# label is struck, and dropped the same way, by an alignment line in `segments` or by the
# axis line of another scale -- scales brought close together put one's labels across the
# other's axis.
function _resolve_tick_collisions!(
    scaledrawings::Vector{ScaleDrawing},
    style::NomogramStyle,
    occupied::Vector{NTuple{4,Float64}},
    segments::Vector{Tuple{Point,Point}},
)
    walls = [
        (j, seg) for (j, sd) in enumerate(scaledrawings) if sd.kind === :scale for
        seg in _segments(sd.path)
    ]
    for level in (1, 2), (j, sd) in enumerate(scaledrawings), i in eachindex(sd.ticks)
        t = sd.ticks[i]
        (t.level == level && !isnothing(t.label)) || continue
        box = label_box(t.labelpos, t.label, style.ticklabelsize, t.align, t.rotation)
        hits(seg) = segment_hits_box(seg[1], seg[2], box, style.collisionpad)
        struck = any(hits, segments) || any(w -> w[1] != j && hits(w[2]), walls)
        if struck || any(b -> _overlaps(box, b, style.collisionpad), occupied)
            sd.ticks[i] = _unlabel(t)
        else
            push!(occupied, box)
        end
    end
    return occupied
end

"""
    drawing(nom::AbstractNomogram; style = NomogramStyle(), isopleths = Isopleth[])

Lay out a nomogram: choose and position ticks, resolve label collisions, and place titles
and any alignment lines. Returns a [`NomogramDrawing`](@ref) of pure geometry, shared by
every rendering backend.
"""
function drawing(
    nom::AbstractNomogram;
    style::NomogramStyle = NomogramStyle(),
    isopleths::AbstractVector{Isopleth} = Isopleth[],
)
    style.ticklabelrotation in (:horizontal, :tangent) || throw(
        ArgumentError(
            "ticklabelrotation must be :horizontal or :tangent, got " *
            ":$(style.ticklabelrotation)",
        ),
    )

    occupied = NTuple{4,Float64}[]
    paths = [_path(s) for s in nom.scales]
    titles = _place_titles(nom, paths, style, occupied)
    scaledrawings = ScaleDrawing[]
    for (s, path, (titlepos, titlealign)) in zip(nom.scales, paths, titles)
        majors, minors = candidate_ticks(s.variable, s.domain)
        ticks = thin_ticks(s, majors, minors, style)
        push!(
            scaledrawings,
            ScaleDrawing(
                name(s),
                axislabel(s.variable),
                path,
                ticks,
                titlepos,
                titlealign,
                s.kind,
            ),
        )
    end

    # Alignment-line geometry is settled before tick labels are, because it does not
    # depend on them: an isopleth's segment comes straight from the solution. Knowing
    # where the lines run lets the label pass reject anything they strike through.
    segments = Tuple{Point,Point}[]
    strikes = Tuple{Point,Point}[]
    for iso in isopleths
        length(iso.points) >= 2 || continue
        a = argmin(p -> p[1], iso.points)
        b = argmax(p -> p[1], iso.points)
        u = _unit(_sub(b, a))
        oh = style.isopleth_overhang
        push!(segments, (_offset(a, u, -oh), _offset(b, u, oh)))
        # Only the span between the outermost crossings counts against labels. The
        # overhang past them is a flourish, and it falls exactly where the outer scales'
        # labels sit -- including the label for the value being read.
        push!(strikes, (a, b))
    end

    if style.resolve_collisions
        _resolve_tick_collisions!(
            scaledrawings,
            style,
            occupied,
            style.avoid_isopleths ? strikes : Tuple{Point,Point}[],
        )
    end

    bb = [Inf, Inf, -Inf, -Inf]
    for sd in scaledrawings
        for p in sd.path
            _expand!(bb, p)
        end
        for t in sd.ticks
            _expand!(bb, _offset(t.position, t.normal, t.length))
            isnothing(t.label) && continue
            _expand_box!(
                bb,
                label_box(t.labelpos, t.label, style.ticklabelsize, t.align, t.rotation),
            )
        end
        isempty(sd.label) && continue
        _expand_box!(
            bb,
            label_box(sd.titlepos, sd.label, style.titlesize, sd.titlealign, 0.0),
        )
    end

    isodrawings = IsoplethDrawing[]
    drawable = [iso for iso in isopleths if length(iso.points) >= 2]
    for (iso, seg) in zip(drawable, segments)
        pts = iso.points
        u = _unit(_sub(seg[2], seg[1]))
        _expand!(bb, seg[1])
        _expand!(bb, seg[2])
        labelpos = _offset(seg[2], (1.0, 0.0), style.ticklabelgap)
        labelalign = (:left, :center)
        if !isempty(iso.label)
            labelpos, box = nudge_clear(
                iso.label,
                labelpos,
                (-u[2], u[1]),
                style.ticklabelsize,
                labelalign,
                style,
                occupied,
            )
            _expand_box!(bb, box)
        end
        push!(
            isodrawings,
            IsoplethDrawing(seg, pts, iso.label, labelpos, labelalign, iso.group),
        )
    end

    titlepos = (nom.width / 2, bb[4] + style.titlegap)
    titlealign = (:center, :bottom)
    if !isempty(nom.title)
        _expand_box!(bb, label_box(titlepos, nom.title, style.titlesize, titlealign, 0.0))
    end

    return NomogramDrawing(
        scaledrawings,
        isodrawings,
        nom.title,
        titlepos,
        titlealign,
        (bb[1], bb[2], bb[3], bb[4]),
        style,
    )
end
