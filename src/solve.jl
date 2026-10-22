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

"""
    solve_on(s::Scale, p, q) -> Float64

The value on scale `s` whose point lies on the line through `p` and `q`.

Works for any scale curve, straight or not, by bracketed root finding over the scale's
domain, so the same solver serves forward and inverse problems and every layout.
"""
function solve_on(s::Scale, p::Point, q::Point)
    if isapprox(p[1], q[1]; atol = 1e-12) && isapprox(p[2], q[2]; atol = 1e-12)
        throw(
            ArgumentError(
                "the two known points coincide at $p, so they do not define a line",
            ),
        )
    end
    lo, hi = s.domain
    g = v -> align_residual(s, v, p, q)
    glo, ghi = g(lo), g(hi)
    iszero(glo) && return lo
    iszero(ghi) && return hi
    if glo * ghi > 0
        throw(
            ArgumentError(
                "no solution for :$(name(s)) within its domain " *
                "($(default_tickformat(lo)), $(default_tickformat(hi))): the " *
                "alignment line misses this scale. The inputs are outside the " *
                "range this chart covers.",
            ),
        )
    end
    return Roots.find_zero(g, (lo, hi), Roots.A42())
end

"""
    evaluate(nom::Nomogram; values...) -> NamedTuple

Solve the relation by alignment. Supply every variable but one and the remaining one is
read off its scale; the result is a `NamedTuple` of all variables.

Because the solver only ever asks "where does this line cross that scale", forward and
inverse problems are the same call — any one of the variables may be the unknown.

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
    given[name(target)] = solve_on(target, p, q)
    return NamedTuple{Tuple(known)}(Tuple(given[n] for n in known))
end

"""
    isopleth(nom::Nomogram; label = "", values...) -> Isopleth

The alignment line for one solution of the relation: solves for the missing variable with
[`evaluate`](@ref), then records where the line meets every scale so it can be drawn.
"""
function isopleth(nom::Nomogram; label::AbstractString = "", values...)
    sol = evaluate(nom; values...)
    pts = [point(s, sol[name(s)]) for s in nom.scales]
    return Isopleth(Dict(pairs(sol)), pts, String(label))
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
