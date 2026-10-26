"""
    Grid(u::Variable, v::Variable, surface)

A two-parameter family of points, `surface(u, v) -> (x, y)`, drawn as a mesh of curves. It
takes the place of one of the three scales an alignment line crosses.

Three collinear points relate three variables. A grid carries two variables in a single
point, so a chart with a grid and two ordinary scales relates four. This is a different
construction from [`compound`](@ref), which chains stages through a shared pivot. A grid can
handle relations that a chain cannot: a pivot works only when the first two variables
combine into one intermediate, while a grid handles a true two-parameter family.

# Deriving a grid

With scale `a` on `x = 0` at height `f(a)` and scale `b` on `x = 1` at height `g(b)`, a line
between them meets `x = x₃` at `(1-x₃)·f(a) + x₃·g(b)`. So the line passes through a grid
point `(x₃, y₃)` exactly when

```
α·f(a) + β·g(b) = γ        with   x₃ = β/(α+β),  y₃ = γ/(α+β)
```

So any relation that can be written this way, with the coefficients `α`, `β`, `γ` depending
on the two grid variables, has a grid. The mesh is truly two-dimensional as long as
`(α : β : γ)` sweeps out a region. If it sweeps out only a curve, the grid collapses to a
line, and the relation needed [`compound`](@ref) instead.

See [`grid_nomogram`](@ref).
"""
struct Grid
    u::Variable
    v::Variable
    surface::Function
end

gridpoint(n::Grid, u, v) = Point(n.surface(float(u), float(v)))

"""
    GridNomogram

A [`Grid`](@ref) and two [`Scale`](@ref)s: four variables related by one alignment line.

Build one with [`grid_nomogram`](@ref). Drawing works as for any other chart: the grid's
strands are scales of `kind = :grid`, labelled with their values, and each family's name is
drawn past the end of its row of labels. How [`evaluate`](@ref) solves the chart depends on
which variable is unknown; see its docstring.
"""
struct GridNomogram <: AbstractNomogram
    grid::Grid
    a::Scale
    b::Scale
    scales::Vector{Scale}
    width::Float64
    height::Float64
    title::String
    relation::Union{Nothing,Function}
end

variablenames(nn::GridNomogram) = [nn.grid.u.name, nn.grid.v.name, name(nn.a), name(nn.b)]

function scale(nn::GridNomogram, n::Symbol)
    name(nn.a) === n && return nn.a
    name(nn.b) === n && return nn.b
    if n === nn.grid.u.name || n === nn.grid.v.name
        throw(
            ArgumentError(
                ":$n is carried by this chart's grid, not by a scale of its own, so it " *
                "has no single curve. The grid's strands are in `scales`.",
            ),
        )
    end
    throw(
        KeyError(
            "no variable :$n in this nomogram (have " *
            join(":" .* string.(variablenames(nn)), ", ") *
            ")",
        ),
    )
end

# Which values of a grid variable get a strand. Explicit ticks win; otherwise the same
# round numbers a scale would be labelled with.
function _strand_values(v::Variable, dom::Tuple{Float64,Float64})
    isnothing(v.ticks) || return sort(collect(v.ticks))
    majors, _ = candidate_ticks(v, dom)
    return sort(majors)
end

# One strand of the mesh: the other grid variable swept with this one held fixed. It is a
# `Scale` of `kind = :grid`, labelled once at its end with the held value and carrying no
# ticks of its own, which is how a grid reads.
function _strand_variable(
    v::Variable,
    held::Float64,
    dom::Tuple{Float64,Float64},
    side::Symbol,
    titleend::Symbol,
)
    label = v.tickformat(held)
    tag = replace(label, '.' => 'p', '-' => 'm')
    return Variable(
        Symbol(v.name, :_, tag),
        LinearScale();
        domain = dom,
        label = label,
        ticks = Float64[],
        side = side,
        titleend = titleend,
    )
end

# Strands of one family usually fan out at one end and crowd together at the other. Label
# them at the end where they are furthest apart: at the crowded end the labels pile up,
# get nudged aside, and end up next to the wrong strand.
function _label_end(held_values, other_dom::Tuple{Float64,Float64}, at)
    spread(t) = begin
        pts = [at(h, t) for h in held_values]
        length(pts) < 2 && return 0.0
        maximum(hypot(p[1] - q[1], p[2] - q[2]) for p in pts for q in pts)
    end
    return spread(other_dom[1]) >= spread(other_dom[2]) ? :start : :stop
end

"""
    grid_nomogram(grid::Grid, a => curve_a, b => curve_b; kwargs...) -> GridNomogram

A four-variable chart: a [`Grid`](@ref) carrying two variables, and two ordinary scales
given as curves, in the manner of [`determinant_nomogram`](@ref).

You can write everything in whatever coordinates suit the geometry, and it is all fitted to
the canvas together. The fit is affine, so it keeps collinear points collinear and cannot
change what the chart computes. The grid's surface is moved by the same fit as the scales.

# Keywords
- `width`, `height`, `title`, `preserve_aspect`, `nsamples`: as for
  [`determinant_nomogram`](@ref).
- `relation`: the algebraic relation `(u, v, a, b) -> residual`, for
  [`max_relation_residual`](@ref) to check the geometry against.

# Example

The general quadratic `a·x² + b·x + c = 0`. Written as `a·x² + b·x = -c`, the coefficients
are `α = x²`, `β = x`, `γ = -c`, so the grid in `(x, c)` is

```
x₃ = x/(x² + x) = 1/(1 + x)        y₃ = -c/(x² + x) = -c/(x·(1 + x))
```

```jldoctest
julia> nom = grid_nomogram(
           Grid(
               Variable(:x; domain = (0.5, 5)),
               Variable(:c; domain = (0, 4)),
               (x, c) -> (1 / (1 + x), -c / (x * (1 + x))),
           ),
           Variable(:a; domain = (0.5, 3)) => a -> (0.0, a),
           Variable(:b; domain = (-6, 0)) => b -> (1.0, b);
           relation = (x, c, a, b) -> a * x^2 + b * x + c,
       );

julia> [round(s.x; digits = 8) for s in evaluate_all(nom; a = 1, b = -3, c = 2)]
2-element Vector{Float64}:
 1.0
 2.0
```
"""
function grid_nomogram(
    grid::Grid,
    spec_a::Pair,
    spec_b::Pair;
    width::Real = 70.0,
    height::Real = 100.0,
    title::AbstractString = "",
    relation = nothing,
    preserve_aspect::Bool = false,
    nsamples::Int = 129,
)
    for v in (grid.u, grid.v)
        isnothing(v.domain) &&
            throw(ArgumentError("grid variable :$(v.name) must declare `domain`"))
    end
    va, fa = first(spec_a), last(spec_a)
    vb, fb = first(spec_b), last(spec_b)
    for v in (va, vb)
        v isa Variable ||
            throw(ArgumentError("each scale must be given as `Variable => curve`"))
        isnothing(v.domain) &&
            throw(ArgumentError("variable :$(v.name) must declare `domain`"))
    end
    du, dv = grid.u.domain, grid.v.domain

    vars = Variable[]
    doms = Tuple{Float64,Float64}[]
    curves = Any[]
    kinds = Symbol[]

    us = _strand_values(grid.u, du)
    vs = _strand_values(grid.v, dv)
    uend = _label_end(us, dv, (h, t) -> gridpoint(grid, h, t))
    vend = _label_end(vs, du, (h, t) -> gridpoint(grid, t, h))
    for held in us
        push!(vars, _strand_variable(grid.u, held, dv, :right, uend))
        push!(doms, dv)
        push!(curves, t -> gridpoint(grid, held, t))
        push!(kinds, :grid)
    end
    for held in vs
        push!(vars, _strand_variable(grid.v, held, du, :left, vend))
        push!(doms, du)
        push!(curves, t -> gridpoint(grid, t, held))
        push!(kinds, :grid)
    end
    for (v, f) in ((va, fa), (vb, fb))
        push!(vars, v)
        push!(doms, v.domain)
        push!(curves, f)
        push!(kinds, :scale)
    end

    placed, fit = _fitted_scales(
        vars,
        doms,
        curves,
        float(width),
        float(height),
        preserve_aspect,
        nsamples,
    )
    scales = Scale[
        Scale(s.variable, s.domain, s.curve, s.side, k) for (s, k) in zip(placed, kinds)
    ]
    sa, sb = scales[end-1], scales[end]
    # The surface has to move with the strands drawn from it, so it takes the same fit.
    placedgrid = Grid(grid.u, grid.v, (u, v) -> fit(gridpoint(grid, u, v)))

    return GridNomogram(
        placedgrid,
        sa,
        sb,
        scales,
        float(width),
        float(height),
        String(title),
        relation,
    )
end

"""
    gridpoint(nn::GridNomogram, u, v) -> Point

Where the pair `(u, v)` falls on the chart's grid.
"""
gridpoint(nn::GridNomogram, u, v) = gridpoint(nn.grid, u, v)

# Twice the signed area of the triangle made by the two scale points and the grid point:
# zero exactly when the three are collinear, as everywhere else in the package.
_grid_residual(nn::GridNomogram, u, v, p::Point, q::Point) =
    _cross(_sub(q, p), _sub(gridpoint(nn, u, v), p))

function _solve_grid_axis(nn::GridNomogram, which::Symbol, held, p::Point, q::Point)
    var = which === :u ? nn.grid.u : nn.grid.v
    lo, hi = var.domain
    g =
        which === :u ? (t -> _grid_residual(nn, t, held, p, q)) :
        (t -> _grid_residual(nn, held, t, p, q))
    a, b = _search_bounds(g, lo, hi)
    roots = Roots.find_zeros(g, a, b)
    isempty(roots) && throw(
        ArgumentError(
            "no solution for :$(var.name) within its domain " *
            "($(default_tickformat(lo)), $(default_tickformat(hi))): the alignment " *
            "line never meets that strand of the grid.",
        ),
    )
    return unique!(sort!(clamp.(roots, lo, hi)))
end

function _grid_solutions(nn::GridNomogram, values)
    names = variablenames(nn)
    un, vn, an, bn = names
    given = Dict{Symbol,Float64}()
    for (k, x) in pairs(values)
        k in names || throw(
            ArgumentError(
                "no variable :$k in this nomogram (have " *
                join(":" .* string.(names), ", ") *
                ")",
            ),
        )
        given[k] = _number(_variable(nn, k), x)
    end
    missingnames = [n for n in names if !haskey(given, n)]
    if length(missingnames) != 1
        throw(
            ArgumentError(
                "supply exactly one fewer value than there are variables: this " *
                "nomogram has 4 and $(length(given)) were given, leaving " *
                "$(length(missingnames)) unknown " *
                "($(join(":" .* string.(missingnames), ", ")))",
            ),
        )
    end
    target = only(missingnames)

    if target === an || target === bn
        # Both grid variables are known, so the grid contributes one fixed point and this
        # is the ordinary three-point problem.
        p = gridpoint(nn, given[un], given[vn])
        other = target === an ? nn.b : nn.a
        q = point(other, given[name(other)])
        s = target === an ? nn.a : nn.b
        roots = solve_all_on(s, p, q)
        isempty(roots) && throw(_no_solution_error(s))
        return [merge(given, Dict(target => r)) for r in roots]
    end

    # One grid variable is unknown: the alignment line is already fixed by the two scales,
    # and the answer is where it crosses the strand for the known grid variable. That
    # strand is a curve, so a line can cross it more than once -- which is exactly what
    # happens on a quadratic chart, where both roots are genuine.
    p = point(nn.a, given[an])
    q = point(nn.b, given[bn])
    which = target === un ? :u : :v
    held = given[target === un ? vn : un]
    roots = _solve_grid_axis(nn, which, held, p, q)
    return [merge(given, Dict(target => r)) for r in roots]
end

_as_namedtuple(nn::GridNomogram, d::Dict{Symbol,Float64}) =
    NamedTuple{Tuple(variablenames(nn))}(Tuple(d[n] for n in variablenames(nn)))

"""
    evaluate(nn::GridNomogram; values...) -> NamedTuple

Solve a grid chart. Supply three of its four variables.

How the chart is solved depends on which variable is unknown. If both grid variables are
known, the grid gives a single point, and the problem is an ordinary three-point alignment.
If one of them is unknown, the two scales already fix the line, and the answer is where that
line crosses the strand for the known grid variable. The strand is a curve, so the line may
cross it more than once. As with other charts, [`evaluate`](@ref) throws an error when there
are several answers, and [`evaluate_all`](@ref) returns them all.
"""
evaluate(nn::GridNomogram; values...) = _with_units(nn, _evaluate(nn; values...))

function _evaluate(nn::GridNomogram; values...)
    sols = _grid_solutions(nn, values)
    if length(sols) > 1
        target = only(n for n in variablenames(nn) if !haskey(values, n))
        throw(
            ArgumentError(
                "the alignment line crosses the grid at $(length(sols)) values of " *
                ":$target ($(join(default_tickformat.(sol[target] for sol in sols), ", ")))" *
                ", so there is no single answer to return. Use `evaluate_all` to get " *
                "them all, or narrow that variable's domain to the branch you want. " *
                "`isopleth` with the same values draws the line with every crossing " *
                "marked.",
            ),
        )
    end
    return _as_namedtuple(nn, only(sols))
end

"""
    evaluate_all(nn::GridNomogram; values...) -> Vector{NamedTuple}

Every solution consistent with the three values given. See [`evaluate`](@ref).
"""
evaluate_all(nn::GridNomogram; values...) =
    [_with_units(nn, sol) for sol in _evaluate_all(nn; values...)]

_evaluate_all(nn::GridNomogram; values...) =
    [_as_namedtuple(nn, d) for d in _grid_solutions(nn, values)]

function _variable(nn::GridNomogram, n::Symbol)
    n === nn.grid.u.name && return nn.grid.u
    n === nn.grid.v.name && return nn.grid.v
    return scale(nn, n).variable
end

"""
    isopleth(nn::GridNomogram, solution::NamedTuple; label = "", group = 0) -> Isopleth
    isopleth(nn::GridNomogram; label = "", group = 0, values...) -> Isopleth

The alignment line for one reading: it meets scale `a`, the grid, and scale `b`. When the
line crosses the known strand more than once, the keyword form returns that one line with
every crossing marked, as for a [`Nomogram`](@ref); see [`isopleth`](@ref).
"""
function isopleth(
    nn::GridNomogram,
    solution::NamedTuple;
    label::AbstractString = "",
    group::Int = 0,
)
    names = variablenames(nn)
    for n in names
        haskey(solution, n) || throw(ArgumentError("solution is missing a value for :$n"))
    end
    solution = _numbers(nn, solution)
    un, vn = names[1], names[2]
    pts = Point[
        point(nn.a, solution[name(nn.a)]),
        gridpoint(nn, solution[un], solution[vn]),
        point(nn.b, solution[name(nn.b)]),
    ]
    vals = Dict{Symbol,Float64}(n => float(solution[n]) for n in names)
    return Isopleth(vals, pts, String(label), group)
end

function isopleth(nn::GridNomogram; label::AbstractString = "", group::Int = 0, values...)
    sols = _evaluate_all(nn; values...)
    isempty(sols) && return isopleth(nn, _evaluate(nn; values...); label, group)
    length(sols) == 1 && return isopleth(nn, only(sols); label, group)
    un, vn = variablenames(nn)[1:2]
    at(sol) = Point[
        point(nn.a, sol[name(nn.a)]),
        gridpoint(nn, sol[un], sol[vn]),
        point(nn.b, sol[name(nn.b)]),
    ]
    return _isopleth_of_all(sols, variablenames(nn), at, label, group)
end

function Base.show(io::IO, nn::GridNomogram)
    print(io, "GridNomogram(grid :", nn.grid.u.name, " x :", nn.grid.v.name)
    print(io, "; scales :", name(nn.a), ", :", name(nn.b), ")")
    return nothing
end

"""
    max_relation_residual(nn::GridNomogram; n = 7) -> Float64

The largest disagreement between a grid chart's geometry and the relation it is meant to
encode.

Sweeps both grid variables over an `n × n` grid and the first scale over `n` values, solves
the second scale geometrically, and checks each answer against the relation. Inputs whose
alignment line misses a scale are skipped; `NaN` comes back when there is no relation to
check or nothing in the sweep could be solved.
"""
function max_relation_residual(nn::GridNomogram; n::Int = 7)
    isnothing(nn.relation) && return NaN
    un, vn, an, bn = variablenames(nn)
    worst = 0.0
    checked = 0
    for u in range(nn.grid.u.domain...; length = n),
        v in range(nn.grid.v.domain...; length = n),
        a in range(nn.a.domain...; length = n)

        sols = nothing
        try
            sols = _evaluate_all(nn; NamedTuple{(un, vn, an)}((u, v, a))...)
        catch err
            err isa ArgumentError || rethrow()
        end
        isnothing(sols) && continue
        for sol in sols
            checked += 1
            worst = max(worst, abs(nn.relation(sol[un], sol[vn], sol[an], sol[bn])))
        end
    end
    return checked == 0 ? NaN : worst
end

# Which family a strand belongs to: strands are named `Symbol(variable, :_, value)`. The
# longer name wins, should one grid variable's name begin with the other's.
function _strand_family(nn::GridNomogram, s::Scale)
    s.kind === :grid || return nothing
    n = string(name(s))
    fams = [f for f in (nn.grid.u, nn.grid.v) if startswith(n, string(f.name, "_"))]
    isempty(fams) && return nothing
    return argmax(f -> length(string(f.name)), fams).name
end

# A family's name goes off the end of its row (or column) of strand labels, the way an
# axis title continues off the end of its axis: the higher end first, then the other.
# With one strand there is no row, and it follows the strand's own label outwards.
function _extra_titles(nn::GridNomogram, titles, paths, style)
    out = Tuple{String,Vector{Any}}[]
    for f in (nn.grid.u, nn.grid.v)
        text = axislabel(f)
        isempty(text) && continue
        idx = [i for (i, s) in enumerate(nn.scales) if _strand_family(nn, s) === f.name]
        isempty(idx) && continue
        boxes = [
            label_box(
                titles[i][1],
                axislabel(nn.scales[i].variable),
                style.titlesize,
                titles[i][2],
                0.0,
            ) for i in idx
        ]
        # With one strand there is no row: follow the strand's own label outwards.
        path = paths[idx[1]]
        c1 = _box_centre(boxes[1])
        away = _unit(_sub(c1, argmin(p -> _dist(p, c1), (path[1], path[end]))))
        push!(out, (text, _row_title_ends(boxes, style, away)))
    end
    return out
end
