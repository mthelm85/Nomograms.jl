# --- concurrency charts: tracing and drawing --------------------------------------------

# Marching squares: the level `c` of the samples `V` on the grid, as polylines through the
# points where it crosses the grid's edges. An edge is keyed `(1, i, j)` from node (i, j)
# to (i+1, j), or `(2, i, j)` from (i, j) to (i, j+1). Each edge is shared by at most two
# cells and each cell links it to one other edge, so the links form simple paths and loops.
#
# `refine(A, B)`, if given, places the crossing on the edge from node A to node B in place of
# linear interpolation -- for a region's edge, where the samples are only in or out.
function _contour(V::Matrix{Float64}, Xs, Ys, c::Float64; refine = nothing)
    nx, ny = size(V)
    links = Dict{NTuple{3,Int},Vector{NTuple{3,Int}}}()
    function link!(a, b)
        push!(get!(links, a, NTuple{3,Int}[]), b)
        push!(get!(links, b, NTuple{3,Int}[]), a)
    end
    for j = 1:(ny-1), i = 1:(nx-1)
        v1, v2, v3, v4 = V[i, j], V[i+1, j], V[i+1, j+1], V[i, j+1]
        (isnan(v1) || isnan(v2) || isnan(v3) || isnan(v4)) && continue
        a1, a2, a3, a4 = v1 > c, v2 > c, v3 > c, v4 > c
        bottom, right, top, left = (1, i, j), (2, i + 1, j), (1, i, j + 1), (2, i, j)
        edges = NTuple{3,Int}[]
        a1 != a2 && push!(edges, bottom)
        a2 != a3 && push!(edges, right)
        a3 != a4 && push!(edges, top)
        a4 != a1 && push!(edges, left)
        if length(edges) == 2
            link!(edges[1], edges[2])
        elseif length(edges) == 4
            # A saddle. If the centre sides with corner 1, corners 1 and 3 join through it
            # and the curve cuts off corners 2 and 4; otherwise it cuts off 1 and 3.
            if ((v1 + v2 + v3 + v4) / 4 > c) == a1
                link!(bottom, right)
                link!(top, left)
            else
                link!(bottom, left)
                link!(right, top)
            end
        end
    end
    function at(e)
        k, i, j = e
        if !isnothing(refine)
            return refine((Xs[i], Ys[j]), k == 1 ? (Xs[i+1], Ys[j]) : (Xs[i], Ys[j+1]))
        end
        a = V[i, j]
        if k == 1
            t = (c - a) / (V[i+1, j] - a)
            return (Xs[i] + t * (Xs[i+1] - Xs[i]), Ys[j])
        else
            t = (c - a) / (V[i, j+1] - a)
            return (Xs[i], Ys[j] + t * (Ys[j+1] - Ys[j]))
        end
    end
    seen = Set{NTuple{3,Int}}()
    function walk(start)
        path = [start]
        push!(seen, start)
        cur = start
        while true
            i = findfirst(e -> !(e in seen), links[cur])
            isnothing(i) && break
            cur = links[cur][i]
            push!(path, cur)
            push!(seen, cur)
        end
        # Back where it started: a closed loop.
        length(path) > 2 && start in links[cur] && push!(path, start)
        return path
    end
    lines = Vector{Point}[]
    keyset = sort!(collect(keys(links)))
    # Open curves first, from one of their ends; whatever is left is a loop.
    for e in keyset
        (e in seen || length(links[e]) != 1) && continue
        push!(lines, at.(walk(e)))
    end
    for e in keyset
        e in seen || push!(lines, at.(walk(e)))
    end
    return lines
end

# The gradient of a page-point function, by central differences, falling back to one side
# where the other has no value.
function _gradient(φ, P::Point; h = 1e-4)
    function d(k)
        e = k == 1 ? (h, 0.0) : (0.0, h)
        fp, fm = φ((P[1] + e[1], P[2] + e[2])), φ((P[1] - e[1], P[2] - e[2]))
        isfinite(fp) && isfinite(fm) && return (fp - fm) / 2h
        f0 = φ(P)
        isfinite(fp) && return (fp - f0) / h
        isfinite(fm) && return (f0 - fm) / h
        return NaN
    end
    return (d(1), d(2))
end

# Move a traced point exactly onto the level `c`: Newton steps along the gradient, or
# along the frame's edge for a point on it, so that a curve still ends exactly at the edge.
# A step that fails, or wanders further than a grid cell, leaves the point where linear
# interpolation put it, which is already close.
function _project(φ, c::Float64, P::Point, W::Float64, H::Float64)
    onx = P[1] <= 0 || P[1] >= W
    ony = P[2] <= 0 || P[2] >= H
    onx && ony && return P
    Q = P
    for _ = 1:3
        v = φ(Q)
        isfinite(v) || return P
        g = _gradient(φ, Q)
        g = onx ? (0.0, g[2]) : ony ? (g[1], 0.0) : g
        n2 = g[1]^2 + g[2]^2
        (isfinite(n2) && n2 > 0) || return P
        Q = (Q[1] - (v - c) * g[1] / n2, Q[2] - (v - c) * g[2] / n2)
    end
    (0 <= Q[1] <= W && 0 <= Q[2] <= H && _dist(Q, P) <= _GRID_STEP) || return P
    return Q
end

# Douglas-Peucker: drop points within `tol` of the line through the points kept around them.
function _simplify(pts::Vector{Point}, tol::Float64)
    length(pts) <= 2 && return pts
    keep = falses(length(pts))
    keep[1] = keep[end] = true
    stack = [(1, length(pts))]
    while !isempty(stack)
        a, b = pop!(stack)
        b - a < 2 && continue
        chord = _sub(pts[b], pts[a])
        len = _norm(chord)
        off(p) = len == 0 ? _dist(p, pts[a]) : abs(_cross(chord, _sub(p, pts[a]))) / len
        k = argmax(i -> off(pts[i]), (a+1):(b-1))
        if off(pts[k]) > tol
            keep[k] = true
            push!(stack, (a, k), (k, b))
        end
    end
    return pts[keep]
end

# Line segments bucketed by the square cells of side `size` they pass through, so that a
# label's box need be tested only against the segments near it. Each line carries an
# owner, and a label is never struck by a line it owns.
struct _SegmentIndex
    segments::Vector{Tuple{Point,Point}}
    owners::Vector{Int}
    cells::Dict{Tuple{Int,Int},Vector{Int}}
    size::Float64
end

function _SegmentIndex(lines; size::Float64 = 4.0)
    segments, owners = Tuple{Point,Point}[], Int[]
    cells = Dict{Tuple{Int,Int},Vector{Int}}()
    cell(x) = floor(Int, x / size)
    for (who, path) in lines, (a, b) in _segments(path)
        push!(segments, (a, b))
        push!(owners, who)
        i = length(segments)
        for cx = cell(min(a[1], b[1])):cell(max(a[1], b[1])),
            cy = cell(min(a[2], b[2])):cell(max(a[2], b[2]))

            push!(get!(cells, (cx, cy), Int[]), i)
        end
    end
    return _SegmentIndex(segments, owners, cells, size)
end

function _struck(ix::_SegmentIndex, box::NTuple{4,Float64}, own::Int)
    cell(x) = floor(Int, x / ix.size)
    for cx = cell(box[1]):cell(box[3]), cy = cell(box[2]):cell(box[4])
        for i in get(ix.cells, (cx, cy), ())
            (own > 0 && ix.owners[i] == own) && continue
            a, b = ix.segments[i]
            segment_hits_box(a, b, box, 0.0) && return true
        end
    end
    return false
end

_length(path::Vector{Point}) =
    sum(_dist(path[k], path[k+1]) for k = 1:(length(path)-1); init = 0.0)

# The values a field gets curves at, each with its level (1 for the roundest, which are
# labelled) and the gap, in value, to its nearest neighbour at least as round as itself.
# A curve at that value stands `gap / |∇φ|` from that neighbour on the page.
function _curve_values(v::Variable)
    majors, minors = candidate_ticks(v, v.domain)
    vals = Tuple{Float64,Int,Tuple{Int,Int}}[]
    for (xs, level) in ((majors, 1), (minors, 2)), x in xs
        any(t -> isapprox(t[1], x; rtol = 1e-9, atol = 1e-12), vals) && continue
        tier, _ = _tick_tier(v, v.domain, float(x), level)
        push!(vals, (float(x), level, (level, tier)))
    end
    out = Tuple{Float64,Int,Float64}[]
    for (x, level, rank) in vals
        others = [y for (y, _, r) in vals if r <= rank && y != x]
        gap = isempty(others) ? Inf : minimum(y -> abs(y - x), others)
        push!(out, (x, level, gap))
    end
    return sort!(out; by = first)
end

# One family's curves on the page: for each value, the pieces of its curve that stand far
# enough from their neighbours. `(value, level, pieces)` per value with anything to draw.
function _trace_field(c::ConcurrencyChart, f::_Field, style::NomogramStyle)
    Xs, Ys = _grid(c)
    V = _samples(c, f, Xs, Ys)
    φ = P -> _field_at(c, f, P)
    W, H = c.width, c.height
    out = Tuple{Float64,Int,Vector{Vector{Point}}}[]
    for (value, level, gap) in _curve_values(f.variable)
        pieces = Vector{Point}[]
        for line in _contour(V, Xs, Ys, value)
            current = Point[]
            for P in line
                P = _project(φ, value, P, W, H)
                room = isfinite(gap) ? gap / hypot(_gradient(φ, P)...) : Inf
                if room >= style.min_curve_spacing
                    push!(current, P)
                else
                    length(current) >= 2 && push!(pieces, current)
                    current = Point[]
                end
            end
            length(current) >= 2 && push!(pieces, current)
        end
        pieces = [_simplify(p, 0.01) for p in pieces if _length(p) >= 1.0]
        isempty(pieces) || push!(out, (value, level, pieces))
    end
    return out
end

# Lines across the chart at an axis's tick values, clipped to the region.
function _gridlines(c::ConcurrencyChart, axis::Symbol, values)
    W, H = c.width, c.height
    span, across = axis === :x ? (W, H) : (H, W)
    topage, _ = _axis_page(axis === :x ? c.x : c.y, span)
    inside(P) = _inside(c, P)
    out = Vector{Point}[]
    for t in values
        p = topage(t)
        1e-9 * span < p < span * (1 - 1e-9) || continue
        n = ceil(Int, across / _GRID_STEP) + 1
        current = Point[]
        for q in range(0.0, across; length = n)
            P = axis === :x ? (p, q) : (q, p)
            if inside(P)
                push!(current, P)
            else
                length(current) >= 2 && push!(out, [current[1], current[end]])
                current = Point[]
            end
        end
        length(current) >= 2 && push!(out, [current[1], current[end]])
    end
    return out
end

# The edge of the chart's region, where it runs inside the frame: traced on the grid and
# placed on each grid edge it crosses by bisection, so it is exact to well within a stroke
# however coarse the grid.
function _region_outline(c::ConcurrencyChart)
    isnothing(c.region) && return Vector{Point}[]
    Xs, Ys = _grid(c)
    V = [_inside(c, (X, Y)) ? 1.0 : 0.0 for X in Xs, Y in Ys]
    function bisect(A, B)
        inA = _inside(c, A)
        for _ = 1:40
            M = ((A[1] + B[1]) / 2, (A[2] + B[2]) / 2)
            _inside(c, M) == inA ? (A = M) : (B = M)
        end
        return ((A[1] + B[1]) / 2, (A[2] + B[2]) / 2)
    end
    lines = _contour(V, Xs, Ys, 0.5; refine = bisect)
    return [_simplify(p, 0.01) for p in lines if length(p) >= 2]
end

# Where a label can go at either end of a traced piece, best first: past the right or top
# edge of the frame, where nothing else is; then continuing the curve from an end inside
# the frame; last past the left or bottom edge, among the axes' own labels.
# `(rank, anchor, align)` for each.
function _label_spots(piece::Vector{Point}, W::Float64, H::Float64, gap::Float64)
    spots = Tuple{Int,Point,NTuple{2,Symbol}}[]
    tol = 1e-6 * max(W, H)
    for (P, Q) in ((piece[end], piece[end-1]), (piece[1], piece[2]))
        if abs(P[1] - W) <= tol
            push!(spots, (1, (W + gap, P[2]), (:left, :center)))
        elseif abs(P[2] - H) <= tol
            push!(spots, (2, (P[1], H + gap), (:center, :bottom)))
        elseif abs(P[1]) <= tol
            push!(spots, (4, (-gap, P[2]), (:right, :center)))
        elseif abs(P[2]) <= tol
            push!(spots, (5, (P[1], -gap), (:center, :top)))
        else
            d = _unit(_sub(P, Q))
            push!(spots, (3, _offset(P, d, gap), align_for(d)))
        end
    end
    return sort!(spots; by = first)
end

"""
    drawing(chart::ConcurrencyChart; style = NomogramStyle(), readings = Reading[])

Lay out a [`ConcurrencyChart`](@ref): its axes and frame, each field's curves and their
labels, and any [`Reading`](@ref)s, each marked with lines dropped to the two axes. Returns
a [`NomogramDrawing`](@ref), which both renderers draw.
"""
function drawing(
    c::ConcurrencyChart;
    style::NomogramStyle = NomogramStyle(),
    readings::AbstractVector{Reading} = Reading[],
    isopleths::AbstractVector{Isopleth} = Isopleth[],
)
    isempty(isopleths) || throw(
        ArgumentError(
            "a concurrency chart is read at a point, not along a line: mark readings " *
            "with `reading` and pass them as `readings = [...]`",
        ),
    )
    W, H = c.width, c.height
    occupied = NTuple{4,Float64}[]
    nowhere = ((0.0, 0.0), (:center, :center))

    # The two axes, with ticks, and the rest of the frame, bare.
    xticks = thin_ticks(c.x, candidate_ticks(c.x.variable, c.x.domain)..., style)
    yticks = thin_ticks(c.y, candidate_ticks(c.y.variable, c.y.domain)..., style)
    frame = [
        ScaleDrawing(name(c.x), "", [(0.0, 0.0), (W, 0.0)], xticks, nowhere..., :scale),
        ScaleDrawing(name(c.y), "", [(0.0, 0.0), (0.0, H)], yticks, nowhere..., :scale),
        ScaleDrawing(:frame_top, "", [(0.0, H), (W, H)], TickMark[], nowhere..., :scale),
        ScaleDrawing(:frame_right, "", [(W, 0.0), (W, H)], TickMark[], nowhere..., :scale),
    ]
    # The region's edge is part of the frame: curves stop at it.
    for (k, path) in enumerate(_region_outline(c))
        push!(
            frame,
            ScaleDrawing(Symbol(:region_, k), "", path, TickMark[], nowhere..., :scale),
        )
    end
    if style.resolve_collisions
        _resolve_tick_collisions!(frame, style, occupied, Tuple{Point,Point}[])
    else
        for sd in frame, t in sd.ticks
            isnothing(t.label) || push!(
                occupied,
                label_box(t.labelpos, t.label, style.ticklabelsize, t.align, 0.0),
            )
        end
    end

    # The axes' titles: centred below the x axis's labels, and above the y axis.
    texts = TextDrawing[]
    walls = [(k, seg) for (k, sd) in enumerate(frame) for seg in _segments(sd.path)]
    lowest = minimum((b[2] for b in occupied); init = -style.majorticklength)
    for (s, ends, own) in (
        (c.x, [((W / 2, lowest - style.ticklabelgap), (:center, :top), (0.0, -1.0))], 1),
        (c.y, [((0.0, H + style.titlegap), (:center, :bottom), (0.0, 1.0))], 2),
    )
        text = axislabel(s.variable)
        isempty(text) && continue
        pos, align =
            style.resolve_collisions ?
            _clear_spot(text, ends, style, occupied, walls, own) : (ends[1][1], ends[1][2])
        push!(texts, TextDrawing(text, pos, align))
    end

    # Lines across at the axes' ticks, then every field's curves.
    curves = ScaleDrawing[]
    if c.grid
        for (axis, ticks) in ((:x, xticks), (:y, yticks))
            lines = _gridlines(c, axis, [t.value for t in ticks if t.level == 1])
            for (k, path) in enumerate(lines)
                push!(
                    curves,
                    ScaleDrawing(
                        Symbol(:grid_, axis, :_, k),
                        "",
                        path,
                        TickMark[],
                        nowhere...,
                        :grid,
                    ),
                )
            end
        end
    end
    traced = [(f, _trace_field(c, f, style)) for f in c.fields]
    # Every line a label must not cross, indexed by where it is, and tagged with the piece
    # it belongs to so a label can sit against its own curve.
    lines = [(0, path) for sd in curves for path in (sd.path,)]
    owner = Dict{Vector{Point},Int}()
    for (_, family) in traced, (_, _, pieces) in family, piece in pieces
        owner[piece] = length(owner) + 1
        push!(lines, (owner[piece], piece))
    end
    for sd in frame
        push!(lines, (0, sd.path))
    end
    index = _SegmentIndex(lines)

    # Labels for the roundest curves, each at the best spot clear of every other label
    # and of every line but its own. One that fits nowhere goes unlabelled.
    labelboxes = Dict{Symbol,Vector{NTuple{4,Float64}}}()
    for (f, family) in traced
        v = f.variable
        boxes = NTuple{4,Float64}[]
        for (value, level, pieces) in family
            text = v.tickformat(value)
            label, labelpos, labelalign, labelled = "", (0.0, 0.0), (:center, :center), 0
            if level == 1
                spots = [
                    (rank, k, anchor, align) for (k, piece) in enumerate(pieces) for
                    (rank, anchor, align) in _label_spots(piece, W, H, style.ticklabelgap)
                ]
                sort!(spots; by = s -> (s[1], -_length(pieces[s[2]])))
                for (_, k, anchor, align) in spots
                    box = label_box(anchor, text, style.titlesize, align, 0.0)
                    clear =
                        !style.resolve_collisions || (
                            !any(b -> _overlaps(box, b, style.collisionpad), occupied) &&
                            !_struck(index, box, owner[pieces[k]])
                        )
                    clear || continue
                    label, labelpos, labelalign, labelled = text, anchor, align, k
                    push!(occupied, box)
                    push!(boxes, box)
                    break
                end
            end
            tag = replace(text, '.' => 'p', '-' => 'm', ' ' => '_', '/' => 'd')
            for (k, piece) in enumerate(pieces)
                mine = k == labelled
                push!(
                    curves,
                    ScaleDrawing(
                        Symbol(v.name, :_, tag, :_, k),
                        mine ? label : "",
                        piece,
                        TickMark[],
                        mine ? labelpos : piece[1],
                        labelalign,
                        :grid,
                    ),
                )
            end
        end
        labelboxes[v.name] = boxes
    end

    # Each family's name, off the end of its row of labels.
    for (f, family) in traced
        text = axislabel(f.variable)
        (isempty(text) || isempty(family)) && continue
        boxes = labelboxes[f.variable.name]
        gap = style.ticklabelgap
        ends = if isempty(boxes)
            # Not a single label fitted: name the family at the end of its longest piece.
            piece = argmax(_length, [p for (_, _, ps) in family for p in ps])
            [(a, al, (1.0, 0.0)) for (_, a, al) in _label_spots(piece, W, H, gap)]
        elseif 2count(b -> b[1] >= W, boxes) >= length(boxes)
            # A column of labels down the right edge: above it, or below it.
            top, bottom = maximum(b -> b[4], boxes), minimum(b -> b[2], boxes)
            [
                ((W + gap, top + style.titlegap / 2), (:left, :bottom), (0.0, 1.0)),
                ((W + gap, bottom - style.titlegap / 2), (:left, :top), (0.0, -1.0)),
            ]
        elseif 2count(b -> b[2] >= H, boxes) >= length(boxes)
            # A row along the top: after it, or before it.
            right, left = maximum(b -> b[3], boxes), minimum(b -> b[1], boxes)
            [
                ((right + 2gap, H + gap), (:left, :bottom), (1.0, 0.0)),
                ((left - 2gap, H + gap), (:right, :bottom), (-1.0, 0.0)),
            ]
        else
            away =
                length(boxes) == 1 ?
                _unit(_sub(_box_centre(boxes[1]), (W / 2, H / 2))) : (1.0, 0.0)
            _row_title_ends(boxes, style, away)
        end
        pos, align =
            style.resolve_collisions ? _clear_spot(text, ends, style, occupied, walls, -1) :
            (ends[1][1], ends[1][2])
        push!(texts, TextDrawing(text, pos, align))
    end

    # Readings: the point, with lines dropped to both axes.
    isodrawings = IsoplethDrawing[]
    for (k, r) in enumerate(readings)
        p = r.point
        # The label goes beside the point, in the first of the four quarters around it, at
        # increasing distances, that is clear of other labels and of the curves.
        labelalign = (:left, :bottom)
        labelpos = _offset(p, (1.0, 1.0), style.ticklabelgap)
        if !isempty(r.label)
            for dist in (1, 2, 4, 8), (dx, dy) in ((1, 1), (-1, 1), (1, -1), (-1, -1))
                at = (
                    p[1] + dx * dist * style.ticklabelgap,
                    p[2] + dy * dist * style.ticklabelgap,
                )
                al = (dx > 0 ? :left : :right, dy > 0 ? :bottom : :top)
                box = label_box(at, r.label, style.ticklabelsize, al, 0.0)
                if !any(b -> _overlaps(box, b, style.collisionpad), occupied) &&
                   !_struck(index, box, -1)
                    labelpos, labelalign = at, al
                    break
                end
            end
            push!(
                occupied,
                label_box(labelpos, r.label, style.ticklabelsize, labelalign, 0.0),
            )
        end
        push!(
            isodrawings,
            IsoplethDrawing(((p[1], 0.0), p), [p], r.label, labelpos, labelalign, k),
        )
        push!(isodrawings, IsoplethDrawing(((0.0, p[2]), p), Point[], "", p, labelalign, k))
    end

    scaledrawings = [frame; curves]
    bb = [0.0, 0.0, W, H]
    for sd in scaledrawings
        foreach(p -> _expand!(bb, p), sd.path)
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
    for t in texts
        _expand_box!(bb, label_box(t.pos, t.text, style.titlesize, t.align, 0.0))
    end
    for iso in isodrawings
        isempty(iso.label) && continue
        _expand_box!(
            bb,
            label_box(iso.labelpos, iso.label, style.ticklabelsize, iso.labelalign, 0.0),
        )
    end
    titlepos = (W / 2, bb[4] + style.titlegap)
    titlealign = (:center, :bottom)
    isempty(c.title) ||
        _expand_box!(bb, label_box(titlepos, c.title, style.titlesize, titlealign, 0.0))

    return NomogramDrawing(
        scaledrawings,
        isodrawings,
        c.title,
        titlepos,
        titlealign,
        (bb[1], bb[2], bb[3], bb[4]),
        style,
        texts,
    )
end

function svg(
    c::ConcurrencyChart;
    style::NomogramStyle = NomogramStyle(),
    readings::AbstractVector{Reading} = Reading[],
    isopleths::AbstractVector{Isopleth} = Isopleth[],
    kwargs...,
)
    return svg(drawing(c; style, readings, isopleths); kwargs...)
end
