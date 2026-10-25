_sub(a::Point, b::Point) = (a[1] - b[1], a[2] - b[2])
_cross(a::Point, b::Point) = a[1] * b[2] - a[2] * b[1]

"""
    align_residual(s::Scale, value, p, q) -> Float64

Signed area (twice the triangle area) of `p`, `q` and `point(s, value)`. Zero exactly when
the three are collinear, so this is the residual whose root [`evaluate`](@ref) seeks.

This is d'Ocagne's determinant criterion written as a cross product. Every layout in this
package is a choice of scale curves that makes this criterion agree with the relation being
charted.
"""
align_residual(s::Scale, value, p::Point, q::Point) =
    _cross(_sub(q, p), _sub(point(s, value), p))

function _finite_at(g, x)
    try
        return isfinite(g(x))
    catch
        return false
    end
end

# A root sitting exactly on the end of a domain -- the ratio at the very corner of an
# N-chart, say -- can land an ulp outside it, because the arithmetic that produced the
# domain and the arithmetic that solves for the root are not the same arithmetic. So search
# a hair wider and clamp back. Only where the curve is still defined, though: a transform
# may have a hard edge exactly at its domain, and stepping over it would throw.
function _search_bounds(g, lo::Float64, hi::Float64)
    pad = 1e-9 * max(hi - lo, abs(lo), abs(hi), 1.0)
    a = _finite_at(g, lo - pad) ? lo - pad : lo
    b = _finite_at(g, hi + pad) ? hi + pad : hi
    return (a, b)
end

function _no_solution_error(s::Scale)
    lo, hi = s.domain
    return ArgumentError(
        "no solution for :$(name(s)) within its domain " *
        "($(default_tickformat(lo)), $(default_tickformat(hi))): the alignment " *
        "line misses this scale. The inputs are outside the range this chart covers.",
    )
end

"""
    solve_all_on(s::Scale, p, q) -> Vector{Float64}

Every value on scale `s` whose point lies on the line through `p` and `q`, in increasing
order. Empty when the line misses the scale.

A straight line can cross a curved scale more than once, so there is not always a single
solution. The usual shortcut, checking only the sign of the residual at the ends of the
domain, would report *no* solution for an even number of crossings and pick one arbitrarily
for an odd number. So this function searches the whole domain instead.

For a straight scale the residual is monotone and this returns at most one value, so charts
with only straight scales never have more than one answer.
"""
function solve_all_on(s::Scale, p::Point, q::Point)
    if isapprox(p[1], q[1]; atol = 1e-12) && isapprox(p[2], q[2]; atol = 1e-12)
        throw(
            ArgumentError(
                "the two known points coincide at $p, so they do not define a line",
            ),
        )
    end
    lo, hi = s.domain
    g = v -> align_residual(s, v, p, q)
    a, b = _search_bounds(g, lo, hi)
    L = hypot(_sub(q, p)...)
    found, res, touch = _crossings(s, g, a, b, L)
    isempty(found) && return Float64[]
    found = sort!([(clamp(v, lo, hi), t) for (v, t) in found])
    merged = _merge_close!(found, s, g, L, res, touch)
    return [t ? _closest_approach(g, v, lo, hi) : v for (v, t) in merged]
end

# Every value in [a, b] where the line crosses the scale, found with a guarantee rather than
# by luck. The residual is a cross product with the fixed vector q - p, so moving a page
# distance d along the curve changes it by at most L*d, L = |q - p|. A stretch whose end
# values add up to more than L times its length therefore cannot hide a crossing and is
# dropped; a sign change brackets one; anything else is halved and looked at again. Two
# crossings far enough apart to be told apart on the page -- more than `res` -- are both
# found. Sampling at fixed points, as a general root finder must, misses a pair of close
# crossings whenever no sample lands in the narrow stretch between them.
function _crossings(s::Scale, g, a::Float64, b::Float64, L::Float64; n0::Int = 64)
    vs = collect(range(a, b; length = n0 + 1))
    ps = [point(s, v) for v in vs]
    gs = [g(v) for v in vs]
    total = sum(hypot(_sub(ps[k+1], ps[k])...) for k = 1:n0)
    res = 1e-7 * max(total, floatmin(Float64))
    touch = 1e-12 * max(total, floatmin(Float64))
    # Each root with whether it was found as a touch rather than a sign change.
    roots = Tuple{Float64,Bool}[(v, false) for (v, gv) in zip(vs, gs) if iszero(gv)]
    stack = [(vs[k], vs[k+1], gs[k], gs[k+1], ps[k], ps[k+1]) for k = 1:n0]
    while !isempty(stack)
        va, vb, ga, gb, pa, pb = pop!(stack)
        if !iszero(ga) && !iszero(gb) && (ga < 0) != (gb < 0)
            push!(roots, (Roots.find_zero(g, (va, vb)), false))
            continue
        end
        vm = (va + vb) / 2
        pm = point(s, vm)
        gm = g(vm)
        iszero(gm) && push!(roots, (vm, false))
        arc = hypot(_sub(pm, pa)...) + hypot(_sub(pb, pm)...)
        # The midpoint estimate of the arc is short of the true length by a little on a
        # curve; the margin covers it.
        abs(ga) + abs(gb) > 1.25 * L * arc && continue
        if arc <= res
            # Close enough to see but no sign change: the line touches the curve here, or
            # misses it by less than anyone could draw. Either way it is one reading.
            m, v = minimum(((abs(ga), va), (abs(gm), vm), (abs(gb), vb)))
            m <= L * touch && push!(roots, (v, true))
            continue
        end
        push!(stack, (va, vm, ga, gm, pa, pm), (vm, vb, gm, gb, pm, pb))
    end
    return roots, res, touch
end

# One reading per crossing. Two candidates are the same crossing found twice -- from both
# sides of a sample, or as touches along a stretch where the line runs alongside the curve
# -- when they are within the search resolution on the page, or the line never parts from
# the curve between them by more than a touch. Two genuine crossings close together part
# by far more than that, and both are kept: evaluate_all answers numerically, not only by
# eye.
function _merge_close!(
    roots::Vector{Tuple{Float64,Bool}},
    s::Scale,
    g,
    L::Float64,
    res::Float64,
    touch::Float64,
)
    kept = Tuple{Float64,Bool}[]
    for (r, t) in roots
        if isempty(kept)
            push!(kept, (r, t))
            continue
        end
        r0, t0 = kept[end]
        near = hypot(_sub(point(s, r), point(s, r0))...) <= res
        alongside = all(u -> abs(g(r0 + u * (r - r0))) <= 2L * touch, (0.25, 0.5, 0.75))
        if near || alongside
            # The same crossing: a sign change pins it exactly, so keep that over a touch.
            t0 && !t && (kept[end] = (r, t))
        else
            push!(kept, (r, t))
        end
    end
    return kept
end

# Where the line comes closest to the curve near a touch: the minimum of |residual|, by
# golden-section search over a bracket grown until both ends are farther away. A touch is
# first found only to within the stretch the line runs alongside the curve, which is far
# wider than the precision a double root deserves.
function _closest_approach(g, v::Float64, lo::Float64, hi::Float64)
    f(x) = abs(g(x))
    h = 1e-9 * max(hi - lo, abs(v), 1.0)
    a, b = max(lo, v - h), min(hi, v + h)
    for _ = 1:60
        (f(a) > f(v) || a == lo) && (f(b) > f(v) || b == hi) && break
        h *= 2
        a, b = max(lo, v - h), min(hi, v + h)
    end
    φ = (sqrt(5) - 1) / 2
    c, d = b - φ * (b - a), a + φ * (b - a)
    fc, fd = f(c), f(d)
    for _ = 1:200
        b - a <= 4eps(max(abs(a), abs(b), 1.0)) && break
        if fc < fd
            b, d, fd = d, c, fc
            c = b - φ * (b - a)
            fc = f(c)
        else
            a, c, fc = c, d, fd
            d = a + φ * (b - a)
            fd = f(d)
        end
    end
    return (a + b) / 2
end

"""
    solve_on(s::Scale, p, q) -> Float64

The single value on scale `s` whose point lies on the line through `p` and `q`.

Throws an error if the line misses the scale, or if it crosses the scale more than once,
which can only happen on a curved scale. Use [`solve_all_on`](@ref) to get every crossing.
"""
function solve_on(s::Scale, p::Point, q::Point)
    roots = solve_all_on(s, p, q)
    isempty(roots) && throw(_no_solution_error(s))
    length(roots) == 1 && return only(roots)
    throw(
        ArgumentError(
            "the alignment line crosses :$(name(s)) at $(length(roots)) values " *
            "($(join(default_tickformat.(roots), ", "))), so there is no single " *
            "answer to return. Use `evaluate_all` / `solve_all_on` to get them all, " *
            "or narrow the scale's domain to the branch you want. `isopleth` with the " *
            "same values draws the line with every crossing marked.",
        ),
    )
end

# Shared front half of `evaluate` and `evaluate_all`: validate the given values, work out
# which variable is unknown, and build the two points that fix the alignment line.
function _alignment(nom::Nomogram, values)
    given = Dict{Symbol,Float64}()
    known = variablenames(nom)
    for (k, v) in pairs(values)
        k in known || throw(
            ArgumentError(
                "no variable :$k in this nomogram (have " *
                join(":" .* string.(known), ", ") *
                ")",
            ),
        )
        given[k] = _number(scale(nom, k).variable, v)
    end
    missingnames = filter(n -> !haskey(given, n), known)
    if length(missingnames) != 1
        throw(
            ArgumentError(
                "supply exactly one fewer value than there are variables: this " *
                "nomogram has $(length(known)) variables and $(length(given)) " *
                "were given, leaving $(length(missingnames)) unknown " *
                "($(join(":" .* string.(missingnames), ", ")))",
            ),
        )
    end
    target = scale(nom, only(missingnames))

    # Two known points define the alignment line. With more than three scales any two
    # will do, and disagreement among the rest would mean the scales are inconsistent.
    anchors = [scale(nom, n) for n in known if n !== name(target)]
    for s in anchors
        v = given[name(s)]
        lo, hi = s.domain
        lo - 1e-9 * max(1.0, abs(lo)) <= v <= hi + 1e-9 * max(1.0, abs(hi)) || throw(
            ArgumentError(
                "value $v for :$(name(s)) is outside its domain " *
                "($(default_tickformat(lo)), $(default_tickformat(hi)))",
            ),
        )
    end
    p = point(anchors[1], given[name(anchors[1])])
    q = point(anchors[2], given[name(anchors[2])])
    return (known, given, target, p, q)
end

_with_root(known, given, target, root) =
    NamedTuple{Tuple(known)}(Tuple(n === name(target) ? root : given[n] for n in known))

# A solution as the caller sees it: each value in its variable's unit, if it has one.
# Everything inside the package works on the plain numbers.
_with_units(nom::AbstractNomogram, sol::NamedTuple) =
    NamedTuple{keys(sol)}(Tuple(_attach(_variable(nom, n).unit, sol[n]) for n in keys(sol)))

# A solution handed back in, possibly with units, as plain numbers.
_numbers(nom::AbstractNomogram, sol::NamedTuple) =
    NamedTuple{keys(sol)}(Tuple(_number(_variable(nom, n), sol[n]) for n in keys(sol)))

# The variable behind a name. A net's own two variables have no scale of their own.
_variable(nom::AbstractNomogram, n::Symbol) = scale(nom, n).variable

"""
    evaluate(nom::Nomogram; values...) -> NamedTuple

Solve the relation by alignment. Supply every variable but one, and the remaining one is
read from its scale. The result is a `NamedTuple` of all the variables.

The solver only ever finds where a line crosses a scale, so forward and inverse problems use
the same call: any one of the variables can be the unknown.

Throws an error if the alignment line crosses the unknown's scale more than once, which can
only happen on a curved scale. Use [`evaluate_all`](@ref) to get every solution.

# Example

```jldoctest
julia> nom = additive_nomogram(
           Variable(:a; domain = (0, 10)),
           Variable(:b; domain = (0, 10)),
           Variable(:s),
       );

julia> round(evaluate(nom; a = 2, s = 9).b; digits = 10)
7.0
```
"""
evaluate(nom::Nomogram; values...) = _with_units(nom, _evaluate(nom; values...))

function _evaluate(nom::Nomogram; values...)
    known, given, target, p, q = _alignment(nom, values)
    return _with_root(known, given, target, solve_on(target, p, q))
end

"""
    evaluate_all(nom::Nomogram; values...) -> Vector{NamedTuple}

Every solution of the relation consistent with the values given, in increasing order of
the unknown.

A straight alignment line can cross a curved scale more than once, and each crossing is a
valid solution. [`evaluate`](@ref) throws an error in that case; `evaluate_all` returns all
of them. For a parallel-layout nomogram the result always has exactly one element.

# Example

```jldoctest
julia> nom = additive_nomogram(
           Variable(:a; domain = (0, 10)),
           Variable(:b; domain = (0, 10)),
           Variable(:s),
       );

julia> length(evaluate_all(nom; a = 2, b = 3))
1
```
"""
evaluate_all(nom::Nomogram; values...) =
    [_with_units(nom, sol) for sol in _evaluate_all(nom; values...)]

function _evaluate_all(nom::Nomogram; values...)
    known, given, target, p, q = _alignment(nom, values)
    roots = solve_all_on(target, p, q)
    isempty(roots) && throw(_no_solution_error(target))
    return [_with_root(known, given, target, r) for r in roots]
end

"""
    isopleth(nom::Nomogram; label = "", values...) -> Isopleth
    isopleth(nom::Nomogram, solution::NamedTuple; label = "") -> Isopleth

The alignment line for one solution of the relation. It records where the line meets every
scale, so that it can be drawn.

The first form solves for the missing variable. If the line crosses that variable's scale
more than once (a curved scale can be crossed twice), every crossing is a solution and all
of them lie on this one line. So the result is that line, with each crossing marked. Its
`values` then hold only the values that fix the line, because the unknown has several
values; get those from [`evaluate_all`](@ref).

The second form takes a solution you already have, to draw one reading on its own:

```julia
isos = [isopleth(nom, sol; label = "z = \$(sol.z)") for sol in evaluate_all(nom; p = -3, q = 2)]
```
"""
function isopleth(
    nom::Nomogram,
    solution::NamedTuple;
    label::AbstractString = "",
    group::Int = 0,
)
    for n in variablenames(nom)
        haskey(solution, n) || throw(
            ArgumentError(
                "solution is missing a value for :$n; it needs one per variable " *
                "(have $(join(":" .* string.(keys(solution)), ", ")))",
            ),
        )
    end
    solution = _numbers(nom, solution)
    pts = [point(s, solution[name(s)]) for s in nom.scales]
    values = Dict{Symbol,Float64}(n => float(solution[n]) for n in variablenames(nom))
    return Isopleth(values, pts, String(label), group)
end

function isopleth(nom::Nomogram; label::AbstractString = "", group::Int = 0, values...)
    sols = _evaluate_all(nom; values...)
    length(sols) == 1 && return isopleth(nom, only(sols); label, group)
    at(sol) = [point(s, sol[name(s)]) for s in nom.scales]
    return _isopleth_of_all(sols, variablenames(nom), at, label, group)
end

# One line marking every crossing. The values that fix the line are shared by all its
# readings; the unknown takes a different value at each crossing, so it is left out.
function _isopleth_of_all(sols, names, at, label::AbstractString, group::Int)
    pts = Point[]
    for sol in sols, p in at(sol)
        seen = any(q -> hypot(p[1] - q[1], p[2] - q[2]) <= 1e-9 * (1 + hypot(q...)), pts)
        seen || push!(pts, p)
    end
    shared(n) = all(sol -> isapprox(sol[n], sols[1][n]; rtol = 1e-12, atol = 1e-12), sols)
    vals = Dict{Symbol,Float64}(n => float(sols[1][n]) for n in names if shared(n))
    return Isopleth(vals, pts, String(label), group)
end

"""
    relation_residual(nom::Nomogram, values::NamedTuple) -> Float64

The residual of the *algebraic* relation the nomogram was built from, independent of the
geometry. Use it to check that a chart encodes the equation it is meant to:
`relation_residual` should be zero for anything [`evaluate`](@ref) returns.

Returns `NaN` for nomograms built directly from scales, which carry no algebraic relation.
"""
function relation_residual(nom::Nomogram, values::NamedTuple)
    isnothing(nom.relation) && return NaN
    values = _numbers(nom, values)
    args = (values[name(s)] for s in nom.scales)
    return float(nom.relation(args...))
end
