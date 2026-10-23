# Rebuild a Variable, keeping every field the caller does not replace.
function rebuild(
    v::Variable;
    transform::ScaleTransform = v.transform,
    domain = v.domain,
    tickstyle::Symbol = v.tickstyle,
)
    return Variable(
        v.name,
        transform;
        domain,
        label = v.label,
        ticks = v.ticks,
        minorticks = v.minorticks,
        tickformat = v.tickformat,
        tickstyle,
        side = v.side,
    )
end

withtransform(v::Variable, t::ScaleTransform) = rebuild(v; transform = t)

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
                "f = [$r0, $r1], but :$(v.name) spans f = [$s0, $s1]. The extent " *
                "of the result axis is fixed by the outer two, so this chart " *
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

"""
    ratio_nomogram(v1, v2, v3; width = 70, height = 100, title = "")

An N-chart (also called a Z-chart) for the ratio relation

```
f1(v1) / f2(v2) = f3(v3)
```

Two vertical scales joined by a diagonal, which is where the result is read. `v3` may
leave `domain = nothing`, in which case the range it can reach is derived from the other
two. All three transforms must be strictly positive over their domains.

# Why not a parallel chart?

[`multiplicative_nomogram`](@ref) already charts ratios, by taking logarithms. The N-chart
earns its place by needing no logarithm at all: the scales stay *linear*, so precision is
even along them and small ratios are not squeezed into a corner. It is the layout to reach
for when the quantities are everyday ones -- a rate as distance over time, a concentration,
a unit conversion -- and the reader wants to interpolate by eye.

# The construction

Put `v1` on `x = 0` at `y1 = a1 * f1`, running up from the bottom-left corner, and `v2` on
`x = W` at `y2 = H - a2 * f2`, running down from the top-right. Take the diagonal between
those two corners, `D(s) = (sW, sH)`. A line through a point on each outer scale meets the
diagonal where

```
s = a1*f1 / (a1*f1 + a2*f2)
```

Writing `r = f1/f2` and `k = a1/a2` that is `s = k*r / (k*r + 1)`, a one-to-one map from
the ratio onto the diagonal: reading `s` back gives `r = s / (k * (1 - s))`.

Note what the two scale factors do *not* have to be. The usual textbook derivation makes
them equal and the ratio falls out directly; here the diagonal's own scale absorbs `k`
instead, so each outer axis is free to be fitted to the canvas on its own. Taking
`a1 = H/max(f1)` and `a2 = H/max(f2)` anchors the diagonal at the two corners, which is
what gives the layout its N.

The diagonal is drawn only over the part its scale actually reaches, `f1/f2` between
`min(f1)/max(f2)` and `max(f1)/min(f2)`, so it stops short of the corners -- those are the
ratios 0 and infinity.

# Example

```jldoctest
julia> nom = ratio_nomogram(
           Variable(:distance; domain = (1, 10)),
           Variable(:time; domain = (1, 10)),
           Variable(:speed),
       );

julia> nom[:speed].domain
(0.1, 10.0)

julia> round(evaluate(nom; distance = 6, time = 3).speed; digits = 10)
2.0
```
"""
function ratio_nomogram(
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
    for (v, lo) in ((v1, p0), (v2, q0))
        lo > 0 || throw(
            ArgumentError(
                "an N-chart divides one scale by the other, so both must be " *
                "strictly positive: :$(v.name) reaches f = $lo on its domain " *
                "$(v.domain). Shift the variable, or use `multiplicative_nomogram` " *
                "if a logarithmic chart suits.",
            ),
        )
    end

    a1 = H / p1
    a2 = H / q1
    k = a1 / a2
    r0, r1 = p0 / q1, p1 / q0
    dom3 = _result_domain(v3, (r0, r1))

    # The diagonal is a projective scale -- s = kr/(kr+1) -- so a ratio spanning decades
    # is spread across it roughly logarithmically. Round numbers picked linearly then sit
    # all in one half: over 0.1 to 10 a linear choice offers nothing between 0.1 and 2,
    # which is most of the axis. Decade-based values fit the geometry, so take them unless
    # the caller has said otherwise.
    if v3.tickstyle === :auto && isnothing(v3.ticks) && r1 / r0 >= 10
        v3 = rebuild(v3; domain = dom3, tickstyle = :log)
    end

    t1, t2, t3 = v1.transform, v2.transform, v3.transform
    c1 = u -> (0.0, a1 * transform(t1, u))
    c2 = v -> (W, H - a2 * transform(t2, v))
    function c3(w)
        kr = k * transform(t3, w)
        s = kr / (kr + 1)
        return (s * W, s * H)
    end

    side1 = v1.side === :auto ? :left : v1.side
    side2 = v2.side === :auto ? :right : v2.side
    side3 = v3.side === :auto ? :left : v3.side

    scales = [
        Scale(v1, v1.domain, c1, side1),
        Scale(v2, v2.domain, c2, side2),
        Scale(v3, dom3, c3, side3),
    ]
    relation = (x, y, z) -> transform(t1, x) / transform(t2, y) - transform(t3, z)
    return Nomogram(scales; width = W, height = H, title, relation)
end
