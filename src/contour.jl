# --- contour charts -------------------------------------------------------------------
#
# A chart read at a point instead of along a line. Two axis variables span a rectangle, and
# every other variable -- a field -- is a function of the two, drawn as its curves of
# constant value. Any two values fix a point where their curves cross; every other value is
# read at that point. Every function of two variables has such a chart, separable or not.
#
# Inside, everything is a function of the page point P = (X, Y): the axes map their values
# linearly in their transforms onto [0, width] and [0, height], and a field is its function
# of the two axis values there. Curves are traced on a grid of samples, moved exactly onto
# their level by Newton steps, and kept only where they stand far enough from their
# neighbours on the page -- which is worked out from the field's gradient, so a family
# thins out where it crowds and runs at full density where it has room.

# One field: its variable and its value as a function of the two axis values.
struct _Field
    variable::Variable
    value::_Formula{2}
end

"""
    ContourChart

A chart read at a point: two axis variables span a rectangle, and every other variable is
drawn as a family of curves of constant value across it. Any two values fix the point where
their curves cross, and every other value is read there. Variables on rings around the
chart are read along a line from its centre through the point.

Build one with [`contour_chart`](@ref). Solve it with [`evaluate`](@ref) and
[`evaluate_all`](@ref), given any two of its variables, and mark a point on it with
[`crossing`](@ref).
"""
struct ContourChart <: AbstractNomogram
    x::Scale
    y::Scale
    fields::Vector{_Field}
    region::Union{Nothing,_Formula{2}}
    grid::Bool
    scales::Vector{Scale}
    width::Float64
    height::Float64
    title::String
    # Page coordinate to axis value, compiled once: tracing calls these millions of times.
    xvalue::FunctionWrapper{Float64,Tuple{Float64}}
    yvalue::FunctionWrapper{Float64,Tuple{Float64}}
    axes::Bool
    rings::Vector{Scale}
    centre::Point
    labelat::Dict{Symbol,Tuple{Symbol,Float64}}
end

"""
    Crossing

A point on a [`ContourChart`](@ref) where curves cross: every variable's value
there, and where it is on the page. Build one with [`crossing`](@ref) and draw it with
`drawing(chart; crossings = [...])`.
"""
struct Crossing
    values::Dict{Symbol,Float64}
    point::Point
    label::String
end

Base.getindex(r::Crossing, n::Symbol) = r.values[n]

const _AxisArg = Union{Variable,Pair{Symbol,<:Tuple}}

# The page coordinate of an axis value, and back: linear in the axis's transform, with the
# domain's ends at 0 and `size`, so the axis always runs low to high left to right (or
# bottom to top) whichever way its transform runs.
function _axis_page(v::Variable, size::Float64)
    t = v.transform
    f0, f1 = transform(t, v.domain[1]), transform(t, v.domain[2])
    topage(x) = size * (transform(t, x) - f0) / (f1 - f0)
    fromp(p) = invtransform(t, f0 + p / size * (f1 - f0))
    return topage, fromp
end

_axis_page(s::Scale, size::Float64) = _axis_page(s.variable, size)

function _axis_value(v::Variable, size::Float64)
    return _axis_value(v.transform, v.domain, size)
end

function _axis_value(t::ScaleTransform, dom, size::Float64)
    f0, f1 = transform(t, dom[1]), transform(t, dom[2])
    return FunctionWrapper{Float64,Tuple{Float64}}(
        p -> Float64(invtransform(t, f0 + p / size * (f1 - f0))),
    )
end

_xvalue(c::ContourChart, X) = c.xvalue(Float64(X))
_yvalue(c::ContourChart, Y) = c.yvalue(Float64(Y))

# Whether a page point is in the chart's region; a region that fails to answer says no.
function _inside(c::ContourChart, P)
    isnothing(c.region) && return true
    return try
        c.region(_xvalue(c, P[1]), _yvalue(c, P[2])) > 0.5
    catch err
        err isa Union{DomainError,ArgumentError,InexactError} || rethrow()
        false
    end
end

# A field's value at a page point, or NaN where it has none: outside the region, or where
# the function fails or returns something that is not a finite number.
function _field_at(c::ContourChart, f::_Field, P)
    x, y = _xvalue(c, P[1]), _yvalue(c, P[2])
    return _guarded(c, f.value, x, y)
end

function _guarded(c::ContourChart, F, x, y)
    v = try
        isnothing(c.region) || c.region(x, y) > 0.5 || return NaN
        F(x, y)
    catch err
        err isa Union{DomainError,ArgumentError,InexactError} || rethrow()
        NaN
    end
    return isfinite(v) ? v : NaN
end

# Every variable's value at a page point, as a function per name.
function _value_at(c::ContourChart, n::Symbol, P)
    n === name(c.x) && return _xvalue(c, P[1])
    n === name(c.y) && return _yvalue(c, P[2])
    _isring(c, n) && return _ring_value(c, _ring(c, n), P)
    return _field_at(c, _field(c, n), P)
end

# A ring is read where a line from the centre through the point meets it, on the point's
# side of the centre: the radial cursor of a Smith chart.
function _ring_value(c::ContourChart, ring::Scale, P)
    d = _sub(P, c.centre)
    _norm(d) > 1e-12 * max(c.width, c.height) || return NaN
    roots = solve_all_on(ring, c.centre, P)
    ahead = [r for r in roots if _dot(_sub(point(ring, r), c.centre), d) > 0]
    return isempty(ahead) ? NaN : first(ahead)
end

_dot(a::Point, b::Point) = a[1] * b[1] + a[2] * b[2]

_ring(c::ContourChart, n::Symbol) = c.rings[findfirst(r -> name(r) === n, c.rings)]
_isring(c::ContourChart, n::Symbol) = any(r -> name(r) === n, c.rings)

function _field(c::ContourChart, n::Symbol)
    i = findfirst(f -> f.variable.name === n, c.fields)
    isnothing(i) && throw(KeyError(n))
    return c.fields[i]
end

variablenames(c::ContourChart) =
    [name(c.x); name(c.y); [f.variable.name for f in c.fields]; [name(r) for r in c.rings]]

function scale(c::ContourChart, n::Symbol)
    n === name(c.x) && return c.x
    n === name(c.y) && return c.y
    _isring(c, n) && return _ring(c, n)
    if any(f -> f.variable.name === n, c.fields)
        throw(
            ArgumentError(
                ":$n is drawn as a family of curves on this chart, not as a scale of its own",
            ),
        )
    end
    throw(
        KeyError(
            "no variable :$n in this chart (have " *
            join(":" .* string.(variablenames(c)), ", ") *
            ")",
        ),
    )
end

function _variable(c::ContourChart, n::Symbol)
    n === name(c.x) && return c.x.variable
    n === name(c.y) && return c.y.variable
    _isring(c, n) && return _ring(c, n).variable
    any(f -> f.variable.name === n, c.fields) && return _field(c, n).variable
    return scale(c, n).variable
end

function Base.show(io::IO, c::ContourChart)
    print(io, "ContourChart(:$(name(c.x)), :$(name(c.y)); curves of ")
    join(io, (":$(f.variable.name)" for f in c.fields), ", ")
    isempty(c.rings) || print(io, "; rings ", join((":$(name(r))" for r in c.rings), ", "))
    isempty(c.title) || print(io, "; \"", c.title, "\"")
    print(io, ")")
end

# --- construction -----------------------------------------------------------------------

"""
    contour_chart(x, y, fields...; kwargs...) -> ContourChart
    contour_chart(F, x, y, [z]; kwargs...) -> ContourChart

A chart read at a point. `x` runs along the bottom and `y` up the left side, each a
[`Variable`](@ref) with a domain or simply `name => (lo, hi)`. Every other variable is a
*field*, given as `variable => (x, y) -> value`, and is drawn as its curves of constant
value. Read the chart by finding where the curves of any two known values cross; every
other value is read at that point.

```julia
moody = contour_chart(
    Variable(:Re, LogScale(); domain = (4e3, 1e8), label = "Reynolds number"),
    Variable(:f, LogScale(); domain = (0.008, 0.1), label = "friction factor"),
    Variable(:roughness; domain = (1e-6, 0.05), label = "ε/D") =>
        (Re, f) -> 3.7 * (10^(-1 / (2sqrt(f))) - 2.51 / (Re * sqrt(f))),
)
evaluate(moody; Re = 1e5, roughness = 1e-4).f
```

A field only needs *some* formula for its value at each point. The Colebrook equation
cannot be solved for the friction factor, but it can be solved for the roughness, so the
chart above is the Moody diagram.

The second form charts a formula `z = F(x, y)` and works for any `F`. When `F` passes the
test [`auto_nomogram`](@ref) applies, the axes take the transforms it recovers, and every
curve of `z` comes out a straight line, all of them parallel; otherwise the axes keep the
variables' own transforms.

Field functions, like transforms, take and return plain numbers, in each variable's unit
if it has one. The formula form calls `F` with quantities when the axes have units, as
[`auto_nomogram`](@ref) does.

# Curves
A field gets curves at round values over the range it takes on the chart, or over its
`domain` if it declares one, chosen as its `ticks` or `tickstyle` say. A field left on
`tickstyle = :auto` gets logarithmic values when its range is positive and spans two
decades or more. Each curve is drawn only where it stands at least `min_curve_spacing`
(see [`NomogramStyle`](@ref)) from its neighbours on the page, judged against the curves
at least as round as it is, so a crowded family thins out locally rather than everywhere.
The roundest curves are labelled where they leave the chart, when there is room.

# Charts on another grid
With `axes = false` the two axis variables are coordinates of a plane rather than scales:
nothing is drawn for them, and the plane keeps its proportions, so a circle stays round.
The Smith chart is drawn this way, on the plane of the reflection coefficient, with
resistance and reactance as fields:

```julia
z(u, v) = (1 + complex(u, v)) / (1 - complex(u, v))
smith = contour_chart(
    :u => (-1, 1),
    :v => (-1, 1),
    Variable(:r; domain = (0, 50), ticks = [0.2, 0.5, 1, 2, 5]) => (u, v) -> real(z(u, v)),
    Variable(:x; domain = (-50, 50), ticks = [-2, -1, 0, 1, 2]) => (u, v) -> imag(z(u, v));
    axes = false,
    region = (u, v) -> u^2 + v^2 <= 1,
    rings = [Variable(:angle; domain = (-180, 180)) => θ -> (1.08cosd(θ), 1.08sind(θ))],
    labelat = (r = :x => 0,),
)
```

# Keywords
- `width = 100`: the width of the chart, in chart units.
- `height`: its height; 75 by default, or with `axes = false`, whatever keeps the plane's
  proportions.
- `title`.
- `region`: `(x, y) -> Bool`, the part of the rectangle the chart covers. Fields have no
  value outside it, so curves stop at its edge and nothing there can be read: the
  saturation line of a psychrometric chart, say. Its edge is drawn.
- `grid = true`: whether to draw lines across the chart at the axes' tick values.
- `axes = true`: whether `x` and `y` are drawn as axes, or are only the plane's coordinates.
- `rings`: scales around the chart, each `variable => value -> (x, y)`, a curve in the
  plane's coordinates. A ring is read along a line from `centre` through the point, and
  takes part in [`evaluate`](@ref) like any other variable.
- `centre`: the point rings are read from, in the plane's coordinates; the middle of the
  plane by default.
- `labelat`: where to label a family inside the chart, as `field = other => value`: each
  curve is labelled where it crosses the curve `other = value`, set along the curve on a
  band of background colour. `labelat = (r = :x => 0,)` labels a Smith chart's resistance
  circles along its real axis. Other families are labelled where their curves end.
"""
function contour_chart(
    x::_AxisArg,
    y::_AxisArg,
    fields::Pair...;
    width::Real = 100.0,
    height::Union{Nothing,Real} = nothing,
    title::AbstractString = "",
    region = nothing,
    grid::Bool = true,
    axes::Bool = true,
    rings = (),
    centre = nothing,
    labelat = (;),
)
    vx, vy = _as_variable(x), _as_variable(y)
    for v in (vx, vy)
        isnothing(v.domain) &&
            throw(ArgumentError("axis variable :$(v.name) must declare a domain"))
    end
    isempty(fields) && throw(ArgumentError("a contour chart needs at least one field"))
    fs = _Field[]
    for p in fields
        v = first(p) isa Symbol ? Variable(first(p)) : first(p)
        v isa Variable || throw(
            ArgumentError(
                "a field is `variable => (x, y) -> value`, got $(repr(first(p)))",
            ),
        )
        push!(fs, _Field(v, _Formula(last(p), 2)))
    end
    ringvars = Variable[_as_variable(first(r)) for r in rings]
    for v in ringvars
        isnothing(v.domain) &&
            throw(ArgumentError("ring variable :$(v.name) must declare a domain"))
    end
    names = [vx.name; vy.name; [f.variable.name for f in fs]; [v.name for v in ringvars]]
    allunique(names) || throw(ArgumentError("variable names must be distinct, got $names"))
    W = float(width)
    span(v) = abs(transform(v.transform, v.domain[2]) - transform(v.transform, v.domain[1]))
    H = float(something(height, axes ? 75.0 : W * span(vy) / span(vx)))
    W > 0 && H > 0 || throw(ArgumentError("width and height must be positive"))
    # Ticks outside the frame: below the x axis and left of the y axis.
    vx = rebuild(vx; side = :left)
    vy = rebuild(vy; side = :left)
    xpage, _ = _axis_page(vx, W)
    ypage, _ = _axis_page(vy, H)
    sx = Scale(vx, vx.domain, v -> (xpage(v), 0.0), :left)
    sy = Scale(vy, vy.domain, v -> (0.0, ypage(v)), :left)
    reg = isnothing(region) ? nothing : _Formula((a, b) -> region(a, b) ? 1.0 : 0.0, 2)
    page(P) = (xpage(P[1]), ypage(P[2]))
    C = page(something(centre, (sum(vx.domain) / 2, sum(vy.domain) / 2)))
    # A ring's ticks and labels go outwards, away from the centre, all the way round.
    rs = Scale[]
    for (v, r) in zip(ringvars, rings)
        f = last(r)
        curve = t -> page(f(t))
        t = sum(v.domain) / 2
        d = _tangent(Scale(v, v.domain, curve, :left), t)
        side = _dot((-d[2], d[1]), _sub(curve(t), C)) > 0 ? :ccw : :cw
        push!(rs, Scale(v, v.domain, curve, side))
    end
    at = Dict{Symbol,Tuple{Symbol,Float64}}()
    for (k, p) in pairs(labelat)
        any(f -> f.variable.name === k, fs) ||
            throw(ArgumentError("labelat names :$k, which is not a field of this chart"))
        first(p) in names || throw(
            ArgumentError("labelat for :$k names :$(first(p)), which is not a variable"),
        )
        at[k] = (first(p), float(last(p)))
    end
    c = ContourChart(
        sx,
        sy,
        fs,
        reg,
        grid,
        [sx; sy; rs],
        W,
        H,
        String(title),
        _axis_value(vx, W),
        _axis_value(vy, H),
        axes,
        rs,
        C,
        at,
    )
    # Settle each field's range and tick style now, so that drawing and solving agree.
    for (i, f) in enumerate(c.fields)
        c.fields[i] = _Field(_resolved_field_variable(c, f), f.value)
    end
    return c
end

function contour_chart(
    F,
    x::_AxisArg,
    y::_AxisArg,
    z::Union{Variable,Symbol} = :z;
    kwargs...,
)
    vx, vy, vz = _as_variable(x), _as_variable(y), _as_variable(z)
    for v in (vx, vy)
        isnothing(v.domain) &&
            throw(ArgumentError("input variable :$(v.name) must declare a domain"))
    end
    F, vz = _with_dimensions(F, Variable[vx, vy], vz)
    Fw = _Formula(F, 2)
    # Where the formula is separable, the transforms that make it so put every curve of z
    # on a straight line. Where it is not, the axes keep their own transforms.
    try
        nom = _auto_two(Fw, vx, vy, vz; verify = false)
        vx = nom[vx.name].variable
        vy = nom[vy.name].variable
        if vz.tickstyle === :auto && isnothing(vz.ticks)
            vz = rebuild(vz; tickstyle = resolved_tickstyle(nom[vz.name].variable))
        end
    catch err
        err isa ArgumentError || rethrow()
    end
    return contour_chart(vx, vy, vz => Fw; kwargs...)
end

# The samples every curve is traced from: a grid over the frame, fine enough that a curve
# turning within a cell is rare, with the field's value at each node.
const _GRID_STEP = 0.5

function _grid(c::ContourChart)
    nx = max(ceil(Int, c.width / _GRID_STEP), 8) + 1
    ny = max(ceil(Int, c.height / _GRID_STEP), 8) + 1
    return collect(range(0.0, c.width; length = nx)),
    collect(range(0.0, c.height; length = ny))
end

_samples(c::ContourChart, f::_Field, Xs, Ys) =
    [_field_at(c, f, (X, Y)) for X in Xs, Y in Ys]

# A field's variable with its range and tick style settled from the values it takes.
function _resolved_field_variable(c::ContourChart, f::_Field)
    v = f.variable
    Xs, Ys = _grid(c)
    vals = filter(isfinite, _samples(c, f, Xs, Ys))
    isempty(vals) && throw(
        ArgumentError(
            ":$(v.name) has no finite value anywhere on the chart; check its formula",
        ),
    )
    lo, hi = extrema(vals)
    lo < hi || throw(ArgumentError(":$(v.name) is constant over the chart"))
    dom = something(v.domain, (lo, hi))
    style = v.tickstyle
    if style === :auto && isnothing(v.ticks) && tickstyle(v.transform) === :linear
        style = dom[1] > 0 && dom[2] / dom[1] >= 100 ? :log : :linear
    end
    return rebuild(v; domain = dom, tickstyle = style)
end

# --- reading -----------------------------------------------------------------------------

# Every point on the chart where the two given values hold, as page points.
function _reading_points(c::ContourChart, given::Dict{Symbol,Float64})
    W, H = c.width, c.height
    order = variablenames(c)
    (a, va), (b, vb) = sort!(collect(given); by = p -> findfirst(==(p[1]), order))
    xn, yn = name(c.x), name(c.y)
    xpage, _ = _axis_page(c.x, W)
    ypage, _ = _axis_page(c.y, H)
    if _isring(c, a) || _isring(c, b)
        # A ring value fixes a line out from the centre; the other value is found along it.
        _isring(c, a) &&
            _isring(c, b) &&
            throw(
                ArgumentError(
                    "two ring values fix only a line through the centre, not a point; give " *
                    "one ring value and one value from the chart itself",
                ),
            )
        (ra, rv), (o, ov) = _isring(c, a) ? ((a, va), (b, vb)) : ((b, vb), (a, va))
        R, C = point(_ring(c, ra), rv), c.centre
        ray(t) = (C[1] + t * (R[1] - C[1]), C[2] + t * (R[2] - C[2]))
        gap(t) = _value_at(c, o, ray(t)) - ov
        return [P for P in (ray(t) for t in _line_roots(gap, 0.0, 1.0)) if _inside(c, P)]
    elseif a === xn && b === yn
        P = (xpage(va), ypage(vb))
        return _inside(c, P) ? [P] : Point[]
    elseif a === xn || a === yn
        # One axis value: where a line across the chart at that value meets the other
        # value's curve.
        along = a === xn ? (t -> (xpage(va), t)) : (t -> (t, ypage(va)))
        f = _field(c, b)
        g(t) = _field_at(c, f, along(t)) - vb
        return [along(t) for t in _line_roots(g, 0.0, a === xn ? H : W)]
    end
    # Two fields: where their curves cross, refined from every cell of a coarse grid over
    # which both change sign.
    fa, fb = _field(c, a), _field(c, b)
    sa = fa.variable.domain[2] - fa.variable.domain[1]
    sb = fb.variable.domain[2] - fb.variable.domain[1]
    r(P) = ((_field_at(c, fa, P) - va) / sa, (_field_at(c, fb, P) - vb) / sb)
    Xs = range(0.0, W; length = 61)
    Ys = range(0.0, H; length = round(Int, 60 * H / W) + 2)
    R = [r((X, Y)) for X in Xs, Y in Ys]
    found = Point[]
    for j = 1:(length(Ys)-1), i = 1:(length(Xs)-1)
        corners = (R[i, j], R[i+1, j], R[i, j+1], R[i+1, j+1])
        any(q -> !all(isfinite, q), corners) && continue
        straddles(k) = minimum(q -> q[k], corners) <= 0 <= maximum(q -> q[k], corners)
        (straddles(1) && straddles(2)) || continue
        P = _newton2(r, ((Xs[i] + Xs[i+1]) / 2, (Ys[j] + Ys[j+1]) / 2), W, H)
        isnothing(P) && continue
        any(Q -> _dist(P, Q) <= 1e-7 * max(W, H), found) || push!(found, P)
    end
    return found
end

# Newton's method in the plane, with the Jacobian by central differences. Returns the
# point it converges to, or nothing; a point that leaves the frame is no reading of it.
function _newton2(r, P::Point, W::Float64, H::Float64)
    tol = 1e-9 * max(W, H)
    h = 1e-6 * max(W, H)
    for _ = 1:40
        f = r(P)
        all(isfinite, f) || return nothing
        fx1, fx0 = r((P[1] + h, P[2])), r((P[1] - h, P[2]))
        fy1, fy0 = r((P[1], P[2] + h)), r((P[1], P[2] - h))
        J = [
            (fx1[1]-fx0[1])/2h (fy1[1]-fy0[1])/2h
            (fx1[2]-fx0[2])/2h (fy1[2]-fy0[2])/2h
        ]
        all(isfinite, J) || return nothing
        abs(LinearAlgebra.det(J)) > 1e-14 || return nothing
        step = J \ [f[1], f[2]]
        P = (P[1] - step[1], P[2] - step[2])
        (-tol <= P[1] <= W + tol && -tol <= P[2] <= H + tol) || return nothing
        hypot(step...) <= 1e-12 * max(W, H) && break
    end
    all(q -> abs(q) <= 1e-9, r(P)) || return nothing
    return (clamp(P[1], 0.0, W), clamp(P[2], 0.0, H))
end

# Roots of g on [lo, hi], from sign changes between samples, each refined by bisection.
function _line_roots(g, lo::Float64, hi::Float64; n::Int = 400)
    ts = range(lo, hi; length = n + 1)
    gs = g.(ts)
    roots = Float64[]
    for k = 1:n
        a, b = gs[k], gs[k+1]
        (isfinite(a) && isfinite(b)) || continue
        if a == 0
            push!(roots, ts[k])
        elseif b != 0 && sign(a) != sign(b)
            push!(roots, Roots.find_zero(g, (ts[k], ts[k+1]), Roots.Bisection()))
        end
    end
    isfinite(gs[end]) && gs[end] == 0 && push!(roots, hi)
    return unique(roots)
end

function _reading_values(c::ContourChart, given, P::Point)
    names = variablenames(c)
    return NamedTuple{Tuple(names)}(
        Tuple(haskey(given, n) ? given[n] : _value_at(c, n, P) for n in names),
    )
end

function _given(c::ContourChart, values)
    known = variablenames(c)
    given = Dict{Symbol,Float64}()
    for (k, v) in pairs(values)
        k in known || throw(
            ArgumentError(
                "no variable :$k in this chart (have " *
                join(":" .* string.(known), ", ") *
                ")",
            ),
        )
        given[k] = _number(_variable(c, k), v)
    end
    length(given) == 2 || throw(
        ArgumentError(
            "a contour chart is read from two values, which fix a point; " *
            "$(length(given)) were given",
        ),
    )
    for (k, v) in given
        (k in (name(c.x), name(c.y)) || _isring(c, k)) || continue
        lo, hi = _variable(c, k).domain
        lo - 1e-9 * max(1.0, abs(lo)) <= v <= hi + 1e-9 * max(1.0, abs(hi)) || throw(
            ArgumentError(
                "value $v for :$k is outside its domain " *
                "($(default_tickformat(lo)), $(default_tickformat(hi)))",
            ),
        )
    end
    return given
end

function _evaluate_all(c::ContourChart; values...)
    given = _given(c, values)
    pts = sort!(_reading_points(c, given))
    if isempty(pts)
        desc = join(("$k = $(default_tickformat(v))" for (k, v) in given), " and ")
        throw(ArgumentError("no point on the chart has $desc"))
    end
    return [_reading_values(c, given, P) for P in pts]
end

"""
    evaluate(chart::ContourChart; values...) -> NamedTuple

Read a [`ContourChart`](@ref): given any two of its variables, find the point where
both hold and return every variable's value there. Throws if there is no such point, or
more than one; [`evaluate_all`](@ref) returns them all.
"""
evaluate(c::ContourChart; values...) = _with_units(c, _evaluate(c; values...))

function _evaluate(c::ContourChart; values...)
    sols = _evaluate_all(c; values...)
    length(sols) == 1 || throw(
        ArgumentError(
            "the two curves cross at $(length(sols)) points on this chart; use " *
            "`evaluate_all` to get every crossing",
        ),
    )
    return only(sols)
end

evaluate_all(c::ContourChart; values...) =
    [_with_units(c, sol) for sol in _evaluate_all(c; values...)]

"""
    crossing(chart::ContourChart; label = "", values...) -> Crossing
    crossing(chart::ContourChart, solution::NamedTuple; label = "") -> Crossing

One point on a [`ContourChart`](@ref), for drawing: `drawing(chart; crossings = [r])`
marks it, with lines dropped to both axes. The first form finds the point from two values,
as [`evaluate`](@ref) does; the second takes a solution you already have, such as one of
several from [`evaluate_all`](@ref).
"""
crossing(c::ContourChart; label::AbstractString = "", values...) =
    _crossing(c, _evaluate(c; values...), label)

crossing(c::ContourChart, sol::NamedTuple; label::AbstractString = "") =
    _crossing(c, _numbers(c, sol), label)

function _crossing(c::ContourChart, sol::NamedTuple, label::AbstractString)
    xpage, _ = _axis_page(c.x, c.width)
    ypage, _ = _axis_page(c.y, c.height)
    P = (xpage(sol[name(c.x)]), ypage(sol[name(c.y)]))
    return Crossing(Dict{Symbol,Float64}(pairs(sol)), P, String(label))
end
