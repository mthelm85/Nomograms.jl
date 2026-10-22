# Rebuild a Variable with a different transform, keeping every other field.
function withtransform(v::Variable, t::ScaleTransform)
    return Variable(
        v.name,
        t;
        domain = v.domain,
        label = v.label,
        ticks = v.ticks,
        minorticks = v.minorticks,
        tickformat = v.tickformat,
        side = v.side,
    )
end

function _result_domain(v::Variable, reachable::Tuple{Float64,Float64})
    r0, r1 = reachable
    if isnothing(v.domain)
        if !hasinverse(v.transform)
            throw(
                ArgumentError(
                    "variable :$(v.name) leaves its domain to be derived, which " *
                    "requires an invertible transform. Either pass " *
                    "`domain = (lo, hi)` or build the transform as " *
                    "`FunctionScale(f, finv)`.",
                ),
            )
        end
        wa = float(invtransform(v.transform, r0))
        wb = float(invtransform(v.transform, r1))
        lo, hi = minmax(wa, wb)
        isfinite(lo) && isfinite(hi) && lo < hi || throw(
            ArgumentError(
                "derived domain for :$(v.name) is degenerate ($lo, $hi); " *
                "pass `domain` explicitly",
            ),
        )
        return (lo, hi)
    end
    s0, s1 = frange(v.transform, v.domain)
    ok =
        isapprox(s0, r0; rtol = 1e-6, atol = 1e-12) &&
        isapprox(s1, r1; rtol = 1e-6, atol = 1e-12)
    if !ok
        suggestion = if hasinverse(v.transform)
            lo, hi = minmax(
                float(invtransform(v.transform, r0)),
                float(invtransform(v.transform, r1)),
            )
            "\nThe consistent domain is ($(default_tickformat(lo)), " *
            "$(default_tickformat(hi)))."
        else
            ""
        end
        throw(
            ArgumentError(
                "domain $(v.domain) of the result variable :$(v.name) is " *
                "inconsistent with the relation: the other two scales reach " *
                "f = [$r0, $r1], but :$(v.name) spans f = [$s0, $s1]. The length " *
                "of the middle axis is fixed by the outer two, so this chart " *
                "cannot be drawn." *
                suggestion *
                "\nTo label only part of the axis, keep the full domain and pass " *
                "`ticks = [...]` instead.",
            ),
        )
    end
    return v.domain
end

"""
    additive_nomogram(v1, v2, v3; width = 70, height = 100, title = "")

A three-parallel-axis nomogram for the separable relation

```
f1(v1) + f2(v2) = f3(v3)
```

where each `fi` is the [`ScaleTransform`](@ref) carried by the corresponding
[`Variable`](@ref). `v3` is the *result* variable and is drawn as the middle axis; it may
leave `domain = nothing`, in which case its range is derived from the other two.

# The construction

Place `v1` on `x = 0` at `y1 = a1 * (f1 - p0)` and `v2` on `x = W` at
`y2 = a2 * (f2 - q0)`, where `[p0, p1]` and `[q0, q1]` are the ranges of `f1` and `f2` over
their domains and `a1 = H/(p1-p0)`, `a2 = H/(q1-q0)` fit each to the canvas height `H`.

A straight line through those two points crosses `x = L*W` at

```
y = (1-L) * y1 + L * y2
```

Choosing `L = a1 / (a1 + a2)` makes `(1-L)*a1 == L*a2 == a1*a2/(a1+a2) =: a3`, so that

```
y = a3 * (f1 + f2) - a3 * (p0 + q0)
```

which is exactly a linear scale for `f3 = f1 + f2` with gradient `a3`. Reading a value off
the middle axis therefore solves the relation. Because `a3 * (1/a1 + 1/a2) == 1`, the
middle axis spans the same `[0, H]` as the outer two: its height and position are
*consequences* of the outer scales, not free parameters.

# Example

```jldoctest
julia> nom = additive_nomogram(
           Variable(:a; domain = (0, 10)),
           Variable(:b; domain = (0, 20)),
           Variable(:total),
       );

julia> nom[:total].domain
(0.0, 30.0)

julia> round(evaluate(nom; a = 3, b = 4).total; digits = 10)
7.0
```
"""
function additive_nomogram(
    v1::Variable,
    v2::Variable,
    v3::Variable;
    width::Real = 70.0,
    height::Real = 100.0,
    title::AbstractString = "",
)
    width > 0 || throw(ArgumentError("width must be positive, got $width"))
    height > 0 || throw(ArgumentError("height must be positive, got $height"))
    for v in (v1, v2)
        isnothing(v.domain) && throw(
            ArgumentError(
                "variable :$(v.name) must declare `domain`; only the result " *
                "variable (the third argument) may leave it to be derived",
            ),
        )
    end
    W, H = float(width), float(height)

    p0, p1 = frange(v1.transform, v1.domain)
    q0, q1 = frange(v2.transform, v2.domain)
    a1 = H / (p1 - p0)
    a2 = H / (q1 - q0)
    lam = a1 / (a1 + a2)
    a3 = a1 * a2 / (a1 + a2)
    r0, r1 = p0 + q0, p1 + q1
    dom3 = _result_domain(v3, (r0, r1))

    t1, t2, t3 = v1.transform, v2.transform, v3.transform
    x3 = lam * W
    c1 = v -> (0.0, a1 * (transform(t1, v) - p0))
    c2 = v -> (W, a2 * (transform(t2, v) - q0))
    c3 = v -> (x3, a3 * (transform(t3, v) - r0))

    side1 = v1.side === :auto ? :left : v1.side
    side2 = v2.side === :auto ? :right : v2.side
    side3 = v3.side === :auto ? (lam < 0.5 ? :right : :left) : v3.side

    scales = [
        Scale(v1, v1.domain, c1, side1),
        Scale(v2, v2.domain, c2, side2),
        Scale(v3, dom3, c3, side3),
    ]
    relation = (x, y, z) -> transform(t1, x) + transform(t2, y) - transform(t3, z)
    return Nomogram(scales; width = W, height = H, title, relation)
end

"""
    multiplicative_nomogram(v1, v2, v3; width = 70, height = 100, title = "")

A three-parallel-axis nomogram for the product relation

```
f1(v1) * f2(v2) = f3(v3)
```

Implemented by wrapping each variable's transform in a logarithm and deferring to
[`additive_nomogram`](@ref), which turns the product into a sum. All three transforms must
be strictly positive over their domains. As with `additive_nomogram`, `v3` is the result
and may leave its domain to be derived.

The logarithm's base cancels out of the geometry, so none is offered.

# Example

```jldoctest
julia> nom = multiplicative_nomogram(
           Variable(:x; domain = (1, 10)),
           Variable(:y; domain = (1, 10)),
           Variable(:z),
       );

julia> nom[:z].domain
(1.0, 100.0)

julia> round(evaluate(nom; x = 2, y = 3).z; digits = 10)
6.0
```
"""
function multiplicative_nomogram(v1::Variable, v2::Variable, v3::Variable; kwargs...)
    L = LogScale(10.0)
    return additive_nomogram(
        withtransform(v1, ComposedScale(L, v1.transform)),
        withtransform(v2, ComposedScale(L, v2.transform)),
        withtransform(v3, ComposedScale(L, v3.transform));
        kwargs...,
    )
end
