_sub(a::Point, b::Point) = (a[1] - b[1], a[2] - b[2])
_cross(a::Point, b::Point) = a[1] * b[2] - a[2] * b[1]

"""
    align_residual(s::Scale, value, p, q) -> Float64

Signed area (twice the triangle area) of `p`, `q` and `point(s, value)`. Zero exactly when
the three are collinear, so this is the residual whose root [`evaluate`](@ref) seeks.

This is d'Ocagne's determinant criterion written as a cross product, and it is the only
thing a nomogram asserts; every layout in this package is just a choice of scale curves
that makes the criterion agree with the relation being charted.
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

A straight line can cross a curved scale more than once, so "the" solution is not always
well defined. Checking only the sign of the residual at the ends of the domain -- the
usual shortcut -- would report *no* solution for an even number of crossings and pick an
arbitrary one of an odd number, so the whole domain is searched instead.

For a straight scale the residual is monotone and this returns at most one value, which is
why parallel-layout nomograms never see the ambiguity.
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
    roots = Roots.find_zeros(g, a, b)
    isempty(roots) && return roots
    return unique!(sort!(clamp.(roots, lo, hi)))
end

"""
    solve_on(s::Scale, p, q) -> Float64

The single value on scale `s` whose point lies on the line through `p` and `q`.

Throws if the line misses the scale, and also if it crosses more than once -- an ambiguity
only curved scales can produce. Use [`solve_all_on`](@ref) to get every crossing.
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
            "or narrow the scale's domain to the branch you want.",
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
        given[k] = float(v)
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

"""
    evaluate(nom::Nomogram; values...) -> NamedTuple

Solve the relation by alignment. Supply every variable but one and the remaining one is
read off its scale; the result is a `NamedTuple` of all variables.

Because the solver only ever asks "where does this line cross that scale", forward and
inverse problems are the same call — any one of the variables may be the unknown.

Throws if the alignment line crosses the unknown's scale more than once, which only a
curved scale can do; [`evaluate_all`](@ref) returns every solution instead.

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
function evaluate(nom::Nomogram; values...)
    known, given, target, p, q = _alignment(nom, values)
    return _with_root(known, given, target, solve_on(target, p, q))
end

"""
    evaluate_all(nom::Nomogram; values...) -> Vector{NamedTuple}

Every solution of the relation consistent with the values given, in increasing order of
the unknown.

A straight alignment line can cross a curved scale more than once, and each crossing is a
genuine solution. [`evaluate`](@ref) refuses to choose between them; this returns the lot.
For a parallel-layout nomogram the result always has exactly one element.

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
function evaluate_all(nom::Nomogram; values...)
    known, given, target, p, q = _alignment(nom, values)
    roots = solve_all_on(target, p, q)
    isempty(roots) && throw(_no_solution_error(target))
    return [_with_root(known, given, target, r) for r in roots]
end

"""
    isopleth(nom::Nomogram; label = "", values...) -> Isopleth
    isopleth(nom::Nomogram, solution::NamedTuple; label = "") -> Isopleth

The alignment line for one solution of the relation: records where the line meets every
scale so it can be drawn.

The first form solves for the missing variable with [`evaluate`](@ref). The second takes a
solution you already have, which is how to draw a chart whose scales a line can cross more
than once:

```julia
isos = [isopleth(nom, sol) for sol in evaluate_all(nom; u = 20, v = 20)]
```
"""
function isopleth(nom::Nomogram, solution::NamedTuple; label::AbstractString = "")
    for n in variablenames(nom)
        haskey(solution, n) || throw(
            ArgumentError(
                "solution is missing a value for :$n; it needs one per variable " *
                "(have $(join(":" .* string.(keys(solution)), ", ")))",
            ),
        )
    end
    pts = [point(s, solution[name(s)]) for s in nom.scales]
    values = Dict{Symbol,Float64}(n => float(solution[n]) for n in variablenames(nom))
    return Isopleth(values, pts, String(label))
end

function isopleth(nom::Nomogram; label::AbstractString = "", values...)
    return isopleth(nom, evaluate(nom; values...); label)
end

"""
    relation_residual(nom::Nomogram, values::NamedTuple) -> Float64

The residual of the *algebraic* relation the nomogram was built from, independent of the
geometry. Useful for checking that a chart really encodes the equation it claims to:
`relation_residual` should vanish on anything [`evaluate`](@ref) returns.

Returns `NaN` for nomograms built directly from scales, which carry no algebraic relation.
"""
function relation_residual(nom::Nomogram, values::NamedTuple)
    isnothing(nom.relation) && return NaN
    args = (values[name(s)] for s in nom.scales)
    return float(nom.relation(args...))
end
