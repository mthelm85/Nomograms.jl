# Rebuild a Variable, keeping every field the caller does not replace.
function rebuild(
    v::Variable;
    transform::ScaleTransform = v.transform,
    domain = v.domain,
    tickstyle::Symbol = v.tickstyle,
    unit = v.unit,
    side::Symbol = v.side,
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
        side,
        titleend = v.titleend,
        maxden = v.maxden,
        unit,
        unitlabel = v.unitlabel,
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
    additive_nomogram(v1, v2, v3; width = 70, height = 100, title = "", max_modulus_ratio = 4)

A parallel-scale nomogram for the separable relation

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
the middle axis therefore solves the relation. Because `a3 * (1/a1 + 1/a2) == 1`, the middle
axis spans the same `[0, H]` as the outer two: its height and position are determined by the
outer scales and are not free parameters.

# Keeping the middle axis clear

The middle axis sits at `L = a1 / (a1 + a2)` of the way across, so it lies close to
whichever outer axis has the smaller gradient, which is the one whose `f` spans more. Over
`[1, 1000]` against `[11, 14]` on logarithmic scales the gradients differ about 28-fold, and
the middle axis lands 3% of the way across, nearly on top of its neighbour. So the steeper
gradient is capped at `max_modulus_ratio` times the other, by shortening that axis and
centring it vertically. An outer axis shifted by `d` shifts the line's crossing of the
middle one by `d` times that axis's share, `L` or `1 - L`, so the middle axis is shifted to
match and every reading is unchanged. The middle axis then lies at least
`1/(1 + max_modulus_ratio)` of the width from either outer axis: a fifth, by default. The
shortened axis only loses resolution it could spare, since it was the one stretched
furthest. Pass `max_modulus_ratio = Inf` for full-height outer axes whatever their ranges.

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
    max_modulus_ratio::Real = 4.0,
)
    width > 0 || throw(ArgumentError("width must be positive, got $width"))
    height > 0 || throw(ArgumentError("height must be positive, got $height"))
    max_modulus_ratio >= 1 ||
        throw(ArgumentError("max_modulus_ratio must be at least 1, got $max_modulus_ratio"))
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
    R = float(max_modulus_ratio)
    a1 = min(H / (p1 - p0), R * H / (q1 - q0))
    a2 = min(H / (q1 - q0), R * H / (p1 - p0))
    # A shortened axis is centred; the middle axis moves by each one's share of the shift.
    d1 = (H - a1 * (p1 - p0)) / 2
    d2 = (H - a2 * (q1 - q0)) / 2
    lam = a1 / (a1 + a2)
    a3 = a1 * a2 / (a1 + a2)
    d3 = (1 - lam) * d1 + lam * d2
    r0, r1 = p0 + q0, p1 + q1
    dom3 = _result_domain(v3, (r0, r1))

    t1, t2, t3 = v1.transform, v2.transform, v3.transform
    x3 = lam * W
    c1 = v -> (0.0, a1 * (transform(t1, v) - p0) + d1)
    c2 = v -> (W, a2 * (transform(t2, v) - q0) + d2)
    c3 = v -> (x3, a3 * (transform(t3, v) - r0) + d3)

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
    multiplicative_nomogram(v1, v2, v3; width = 70, height = 100, title = "", max_modulus_ratio = 4)

A parallel-scale nomogram, on logarithmic scales, for the product relation

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

An N chart (or Z chart, by its shape) for the ratio relation

```
f1(v1) / f2(v2) = f3(v3)
```

It has two vertical scales joined by a diagonal, on which the result is read. `v3` may leave
`domain = nothing`, in which case the range it can reach is derived from the other two. All
three transforms must be strictly positive over their domains.

# Compared with a parallel chart

[`multiplicative_nomogram`](@ref) can also chart ratios, by taking logarithms. The N-chart
needs no logarithms: the scales stay *linear*, so precision is even along them and small
ratios are not squeezed into a corner. Use it for everyday quantities, such as a rate as
distance over time, a concentration or a unit conversion, when the reader will want to
interpolate by eye.

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

The two scale factors do not have to be equal. The usual textbook derivation makes them
equal so that the ratio comes out directly. Here the diagonal's own scale absorbs `k`
instead, so each outer axis can be fitted to the canvas separately. Taking `a1 = H/max(f1)`
and `a2 = H/max(f2)` anchors the diagonal at the two corners, which gives the layout its N
shape.

The diagonal is drawn only over the part its scale reaches, `f1/f2` between
`min(f1)/max(f2)` and `max(f1)/min(f2)`, so it stops short of the corners, which correspond
to the ratios 0 and infinity.

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
    max_modulus_ratio::Real = 4.0,
)
    width > 0 || throw(ArgumentError("width must be positive, got $width"))
    height > 0 || throw(ArgumentError("height must be positive, got $height"))
    max_modulus_ratio >= 1 ||
        throw(ArgumentError("max_modulus_ratio must be at least 1, got $max_modulus_ratio"))
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

"""
    concurrent_nomogram(v1, v2, v3; width = 100, height = 50, title = "")

A concurrent-scale nomogram for the relation

```
1/f1(v1) + 1/f2(v2) = 1/f3(v3)
```

which covers resistors in parallel, capacitors in series and the thin-lens equation. Its
three scales start from one point, the zero of every scale: `v1` and `v2` fan out towards
the top corners, and `v3` runs up the middle between them. `v3` may leave
`domain = nothing`, in which case the range it can reach is derived from the other two. All
three transforms must be zero or more over their domains.

# The construction

Put the common zero at `(W/2, 0)`, `v1` at `(W/2 - c·f1, d·f1)`, `v2` at
`(W/2 + c·f2, d·f2)` and `v3` at `(W/2, 2d·f3)`. A line through the first two crosses
`x = W/2` at height `2d·f1·f2 / (f1 + f2)`, which is `2d·f3` exactly when the relation
holds. `c` and `d` are chosen so that the longer outer scale ends at a top corner.

Doerfler derives the same chart with the outer scales at an angle `A` either side of the
middle one and scale factors `m1 = m2 = m3 / (2 cos A)` (see the Bibliography). Stretching
a chart to fill its canvas cannot change what it computes, and it absorbs the angle, so
here `width` and `height` set the spread instead.

# Example

```jldoctest
julia> nom = concurrent_nomogram(
           Variable(:R1; domain = (0, 1000)),
           Variable(:R2; domain = (0, 1000)),
           Variable(:R),
       );

julia> round(evaluate(nom; R1 = 400, R2 = 600).R; digits = 10)
240.0
```
"""
function concurrent_nomogram(
    v1::Variable,
    v2::Variable,
    v3::Variable;
    width::Real = 100.0,
    height::Real = 50.0,
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
    t1, t2, t3 = v1.transform, v2.transform, v3.transform
    p0, p1 = frange(t1, v1.domain)
    q0, q1 = frange(t2, v2.domain)
    for (v, lo) in ((v1, p0), (v2, q0))
        lo >= 0 || throw(
            ArgumentError(
                "the scales of a concurrent-scale nomogram start from a common zero, so " *
                "their values must be zero or more: :$(v.name) reaches f = $lo on its " *
                "domain $(v.domain)",
            ),
        )
    end
    F = max(p1, q1)
    F > 0 || throw(ArgumentError("both outer scales are zero over their whole domains"))
    # f1·f2/(f1 + f2) rises with both, so the corners of the input ranges bound it.
    harmonic(a, b) = a + b == 0 ? 0.0 : a * b / (a + b)
    dom3 = _result_domain(v3, (harmonic(p0, q0), harmonic(p1, q1)))
    c, d = (W / 2) / F, H / F

    scales = [
        Scale(
            v1,
            v1.domain,
            u -> (W / 2 - c * transform(t1, u), d * transform(t1, u)),
            v1.side === :auto ? :left : v1.side,
        ),
        Scale(
            v2,
            v2.domain,
            v -> (W / 2 + c * transform(t2, v), d * transform(t2, v)),
            v2.side === :auto ? :right : v2.side,
        ),
        Scale(
            v3,
            dom3,
            w -> (W / 2, 2d * transform(t3, w)),
            v3.side === :auto ? :right : v3.side,
        ),
    ]
    # In this form the relation has no infinities where a value is zero.
    relation = function (x, y, z)
        a, b, r = transform(t1, x), transform(t2, y), transform(t3, z)
        return r * (a + b) - a * b
    end
    return Nomogram(scales; width = W, height = H, title, relation)
end

"""
    proportional_nomogram(v1, v2, v3, v4; width = 70, height = 100, title = "")

A proportional nomogram for the four-variable relation

```
f1(v1) / f2(v2) = f3(v3) / f4(v4)
```

such as the ideal gas law written `V/T = nR/p`, or the law of sines. It is two N charts
sharing their diagonal. `v1` and `v3` share the left-hand line, with `v1`'s labels outside
and `v3`'s inside; `v2` and `v4` share the right-hand one. Lay a line from `v1` to `v2`,
note where it crosses the diagonal, then lay a second line through that point from `v3`
to `v4`. Any three of the four values give the fourth. All four transforms must be strictly
positive over their domains, and every variable must declare `domain`.

The diagonal carries the common ratio and is drawn bare: it is a pivot line, as in any
[`CompoundNomogram`](@ref), which is what this returns. Following Doerfler (see the
Bibliography), the scale factors keep `m1/m2 = m3/m4`, so the two N charts read the ratio
at the same point of the diagonal.

# Example

```jldoctest
julia> nom = proportional_nomogram(
           Variable(:V; domain = (1, 10)),
           Variable(:T; domain = (250, 400)),
           Variable(:nR; domain = (0.01, 0.2)),
           Variable(:p; domain = (0.5, 5)),
       );

julia> s = evaluate(nom; V = 5, T = 300, p = 2);

julia> round(s.nR; digits = 10)
0.0333333333
```
"""
function proportional_nomogram(
    v1::Variable,
    v2::Variable,
    v3::Variable,
    v4::Variable;
    width::Real = 70.0,
    height::Real = 100.0,
    title::AbstractString = "",
)
    width > 0 || throw(ArgumentError("width must be positive, got $width"))
    height > 0 || throw(ArgumentError("height must be positive, got $height"))
    vs = (v1, v2, v3, v4)
    for v in vs
        isnothing(v.domain) && throw(
            ArgumentError(
                "variable :$(v.name) must declare `domain`; a proportional nomogram has " *
                "no result scale to derive one for",
            ),
        )
    end
    W, H = float(width), float(height)
    ranges = [frange(v.transform, v.domain) for v in vs]
    for (v, (lo, _)) in zip(vs, ranges)
        lo > 0 || throw(
            ArgumentError(
                "a proportional nomogram divides one scale by another, so all four must " *
                "be strictly positive: :$(v.name) reaches f = $lo on its domain $(v.domain)",
            ),
        )
    end
    (p0, p1), (q0, q1), (r0, r1), (s0, s1) = ranges
    # The first pair fill the height, as in an N chart. The second pair keep the same
    # ratio of scale factors, so both read the diagonal alike, and are as large as fits.
    a1, a2 = H / p1, H / q1
    k = a1 / a2
    a3 = min(H / r1, k * H / s1)
    a4 = a3 / k
    lo = max(p0 / q1, r0 / s1)
    hi = min(p1 / q0, r1 / s0)
    lo < hi || throw(
        ArgumentError(
            "the two ratios share no values: f1/f2 runs over " *
            "($(default_tickformat(p0 / q1)), $(default_tickformat(p1 / q0))) and f3/f4 " *
            "over ($(default_tickformat(r0 / s1)), $(default_tickformat(r1 / s0)))",
        ),
    )
    pivot = Variable(:_ratio; domain = (lo, hi), label = "", ticks = Float64[])
    diag = r -> (s = k * r / (k * r + 1); (s * W, s * H))
    t1, t2, t3, t4 = (v.transform for v in vs)
    side(v, default) = v.side === :auto ? default : v.side
    stage1 = Nomogram(
        [
            Scale(v1, v1.domain, u -> (0.0, a1 * transform(t1, u)), side(v1, :left)),
            Scale(v2, v2.domain, u -> (W, H - a2 * transform(t2, u)), side(v2, :right)),
            Scale(pivot, pivot.domain, diag, :left),
        ];
        width = W,
        height = H,
        relation = (x, y, r) -> transform(t1, x) / transform(t2, y) - r,
    )
    stage2 = Nomogram(
        [
            Scale(v3, v3.domain, u -> (0.0, a3 * transform(t3, u)), side(v3, :right)),
            Scale(v4, v4.domain, u -> (W, H - a4 * transform(t4, u)), side(v4, :left)),
            Scale(pivot, pivot.domain, diag, :left),
        ];
        width = W,
        height = H,
        relation = (x, y, r) -> transform(t3, x) / transform(t4, y) - r,
    )
    return compound(stage1, stage2; title)
end

_curvepoint(f, t) = (p = f(t); (float(p[1]), float(p[2])))

# Which side of a scale faces away from the middle of the chart. With parallel axes the
# answer is obvious and the constructions hard-code it; for curves in general position it
# is worth working out, so labels radiate outwards instead of piling into the centre.
function _outward_side(f, domain::Tuple{Float64,Float64}, centre::Point)
    lo, hi = domain
    mid = (lo + hi) / 2
    h = max(1e-6 * (hi - lo), 1e-9)
    a = _curvepoint(f, max(mid - h, lo))
    b = _curvepoint(f, min(mid + h, hi))
    t = (b[1] - a[1], b[2] - a[2])
    n = hypot(t[1], t[2])
    n == 0 && return :right
    left, right = side_normals((t[1] / n, t[2] / n))
    p = _curvepoint(f, mid)
    away = (p[1] - centre[1], p[2] - centre[2])
    # `side` is recorded as :left/:right and resolved against the tangent at draw time, so
    # name whichever normal points away from the middle of the chart.
    return (right[1] * away[1] + right[2] * away[2]) >=
           (left[1] * away[1] + left[2] * away[2]) ? :right : :left
end

"""
    determinant_nomogram(var => curve, var => curve, var => curve, ...; kwargs...)

A nomogram from scales given directly as curves, fitted to the canvas for you.

This is the general constructor. Every three-variable nomogram is a set of parametric curves
positioned so that d'Ocagne's determinant

```
| x₁(u)  y₁(u)  1 |
| x₂(v)  y₂(v)  1 |  =  0
| x₃(w)  y₃(w)  1 |
```

vanishes exactly when the relation holds. [`additive_nomogram`](@ref) and
[`ratio_nomogram`](@ref) are two particular solutions of that condition. This constructor
accepts any other solution you can write down, including curved and concurrent layouts that
parallel axes cannot express.

Each argument is a [`Variable`](@ref) paired with a function from its value to a point. The
curves can use whatever coordinates suit the geometry: their combined bounding box is
measured and mapped onto `width` × `height`. An affine map keeps collinear points collinear,
so the fit cannot change what the chart computes. This includes the default fit, which
scales the two axes by different factors to fill the canvas.

# Keywords
- `width`, `height`: the canvas to fit to.
- `preserve_aspect`: fit with one scale factor and centre the result, rather than
  stretching each axis to fill.
- `relation`: the algebraic relation, if you have it, for [`relation_residual`](@ref) and
  [`max_relation_residual`](@ref) to check the geometry against. It is optional; the chart
  draws and solves without it.
- `title`, `nsamples`: chart title, and how finely curves are sampled when measuring the
  bounding box.

Every variable must have a `domain`. There is no designated result variable, so no domain
can be derived.

# Example

The parallel-axis construction, written out by hand in its natural coordinates. Scales at
`x = 0`, `2` and `1`, with the middle one at half scale, encode `u + v = w`; the fit then
stretches that `2 × 10` sketch onto the full canvas.

```jldoctest
julia> nom = determinant_nomogram(
           Variable(:u; domain = (0, 10)) => u -> (0.0, u),
           Variable(:v; domain = (0, 10)) => v -> (2.0, v),
           Variable(:w; domain = (0, 20)) => w -> (1.0, w / 2);
           relation = (u, v, w) -> u + v - w,
       );

julia> round(evaluate(nom; u = 3, v = 4).w; digits = 10)
7.0

julia> point(nom[:v], 10)
(70.0, 100.0)
```
"""
function determinant_nomogram(
    specs::Pair...;
    width::Real = 70.0,
    height::Real = 100.0,
    title::AbstractString = "",
    relation = nothing,
    preserve_aspect::Bool = false,
    nsamples::Int = 129,
)
    length(specs) >= 3 ||
        throw(ArgumentError("a nomogram needs at least 3 scales, got $(length(specs))"))
    width > 0 || throw(ArgumentError("width must be positive, got $width"))
    height > 0 || throw(ArgumentError("height must be positive, got $height"))
    W, H = float(width), float(height)

    vars = Variable[]
    curves = Any[]
    for spec in specs
        v, f = first(spec), last(spec)
        v isa Variable || throw(
            ArgumentError(
                "each argument must be a `Variable => curve` pair, got a " *
                "$(typeof(v)) on the left",
            ),
        )
        isnothing(v.domain) && throw(
            ArgumentError(
                "variable :$(v.name) must declare `domain`: a determinant nomogram " *
                "has no distinguished result to derive one from",
            ),
        )
        push!(vars, v)
        push!(curves, f)
    end

    doms = Tuple{Float64,Float64}[v.domain for v in vars]
    scales, _ = _fitted_scales(vars, doms, curves, W, H, preserve_aspect, nsamples)
    return Nomogram(scales; width = W, height = H, title, relation)
end

# Measure the curves, map their combined bounding box onto the canvas, and build the
# scales. Returns the fit alongside them, because a caller may hold geometry that is not a
# scale -- a Grid's surface -- and has to move it the same way. Shared by
# `determinant_nomogram`, `project` and `grid_nomogram`, which all start from curves in
# coordinates of their own choosing.
function _fitted_scales(
    vars::Vector{Variable},
    domains::Vector{Tuple{Float64,Float64}},
    curves::Vector,
    W::Float64,
    H::Float64,
    preserve_aspect::Bool,
    nsamples::Int,
)
    # The domain is passed rather than read off the variable: a derived result variable
    # leaves `domain` unset, and it is its scale that knows the range.
    sampled = Vector{Vector{Point}}(undef, length(vars))
    for (i, (v, dom, f)) in enumerate(zip(vars, domains, curves))
        lo, hi = dom
        pts = [_curvepoint(f, t) for t in range(lo, hi; length = nsamples)]
        all(p -> isfinite(p[1]) && isfinite(p[2]), pts) || throw(
            ArgumentError("the curve for :$(v.name) is not finite everywhere on $dom"),
        )
        sampled[i] = pts
    end

    xs = (p[1] for pts in sampled for p in pts)
    ys = (p[2] for pts in sampled for p in pts)
    x0, x1 = extrema(xs)
    y0, y1 = extrema(ys)
    for (span, axis) in ((x1 - x0, "horizontally"), (y1 - y0, "vertically"))
        span > 0 || throw(
            ArgumentError(
                "the scales have no extent $axis, so there is nothing to fit to the " *
                "canvas; check the curves are not all on one line",
            ),
        )
    end

    ax, ay = W / (x1 - x0), H / (y1 - y0)
    if preserve_aspect
        ax = ay = min(ax, ay)
    end
    bx = (W - ax * (x1 - x0)) / 2 - ax * x0
    by = (H - ay * (y1 - y0)) / 2 - ay * y0
    fit(p) = (ax * p[1] + bx, ay * p[2] + by)

    centre = (W / 2, H / 2)
    scales = Scale[]
    for (v, dom, f) in zip(vars, domains, curves)
        placed = t -> fit(_curvepoint(f, t))
        side = v.side === :auto ? _outward_side(placed, dom, centre) : v.side
        push!(scales, Scale(v, dom, placed, side))
    end
    return (scales, fit)
end

"""
    max_relation_residual(nom::Nomogram; n = 9) -> Float64

The largest disagreement between a nomogram's geometry and the relation it is meant to
encode, over a grid of `n × n` inputs.

Getting the determinant form right is the hardest part of designing a nomogram by hand, and
a chart with a small error in it still draws normally. This function solves for the third
variable *geometrically* from the first two, then uses [`relation_residual`](@ref) to check
whether the answer satisfies the relation. A result far from zero means the curves and the
formula disagree.

Returns `NaN` for a nomogram built without a `relation`. Inputs whose alignment line misses
the third scale are skipped, since they say nothing about correctness.
"""
function max_relation_residual(nom::Nomogram; n::Int = 9)
    isnothing(nom.relation) && return NaN
    a, b = nom.scales[1], nom.scales[2]
    target = name(nom.scales[3])
    worst = 0.0
    checked = 0
    for u in range(a.domain...; length = n), v in range(b.domain...; length = n)
        sols = nothing
        try
            sols = _evaluate_all(nom; NamedTuple{(name(a), name(b))}((u, v))...)
        catch err
            # An alignment line that misses the third scale says nothing about whether
            # the formula and the curves agree.
            err isa ArgumentError || rethrow()
        end
        isnothing(sols) && continue
        # Every crossing is a solution and every one has to satisfy the relation, which
        # is the whole point on a chart whose third scale a line can meet twice.
        for sol in sols
            checked += 1
            worst = max(worst, abs(relation_residual(nom, sol)))
        end
    end
    # Vacuous success is worse than no answer: say so rather than returning a clean zero.
    return checked == 0 ? NaN : worst
end

"""
    projective_transform(src, dst) -> Matrix{Float64}

The 3×3 homography taking the four points `src` to the four points `dst`, each given as a
4-tuple (or vector) of `(x, y)`.

Four correspondences fix a projective map exactly, so this is a practical way to describe
one: you say where four corners should land instead of writing out nine coefficients. No
three of the four source points may be collinear, and likewise for the destination points.
"""
function projective_transform(src, dst)
    length(src) == 4 && length(dst) == 4 ||
        throw(ArgumentError("need exactly 4 source and 4 destination points"))
    A = zeros(Float64, 8, 8)
    b = zeros(Float64, 8)
    for i = 1:4
        x, y = float(src[i][1]), float(src[i][2])
        X, Y = float(dst[i][1]), float(dst[i][2])
        A[2i-1, :] .= (x, y, 1, 0, 0, 0, -x * X, -y * X)
        A[2i, :] .= (0, 0, 0, x, y, 1, -x * Y, -y * Y)
        b[2i-1] = X
        b[2i] = Y
    end
    h = try
        A \ b
    catch err
        err isa LinearAlgebra.SingularException || rethrow()
        throw(
            ArgumentError(
                "those four points do not fix a projective map; check that no three " *
                "of the source points are collinear, nor three of the destinations",
            ),
        )
    end
    H = [h[1] h[2] h[3]; h[4] h[5] h[6]; h[7] h[8] 1.0]
    all(isfinite, H) || throw(
        ArgumentError(
            "those four points do not fix a projective map; check that no three of " *
            "the source points are collinear, nor three of the destinations",
        ),
    )
    return H
end

function _project_point(H::AbstractMatrix, p, name::Symbol)
    x, y = p
    w = H[3, 1] * x + H[3, 2] * y + H[3, 3]
    if abs(w) < 1e-12
        throw(
            ArgumentError(
                "the projective map sends a point of :$name to infinity: its vanishing " *
                "line crosses that scale's domain. Move the line clear of the chart, " *
                "or narrow the domain.",
            ),
        )
    end
    return (
        (H[1, 1] * x + H[1, 2] * y + H[1, 3]) / w,
        (H[2, 1] * x + H[2, 2] * y + H[2, 3]) / w,
    )
end

"""
    project(nom::Nomogram, H::AbstractMatrix; kwargs...) -> Nomogram
    project(nom::Nomogram, src => dst; kwargs...) -> Nomogram

Reshape a nomogram with a projective map. The reshaped chart gives exactly the same answers
as the original.

A projective map sends straight lines to straight lines, so it preserves the collinearity
that a nomogram depends on. Everything else about the picture can change. For example, you
can pull open a region where the ticks are crowded together, and the chart still computes
exactly what it did before.

The second form builds the map from four point correspondences with
[`projective_transform`](@ref), which is usually the easiest way to specify it: give the
corners of a quadrilateral in the current chart and where they should end up.

A projective map can move the bounding box anywhere, so the result is fitted to the canvas
again. That fit is affine, and an affine map composed with a projective one is still
projective, so the answers are still unchanged.

# Keywords
- `width`, `height`: canvas for the reshaped chart; defaults to the original's.
- `preserve_aspect`, `nsamples`: as for [`determinant_nomogram`](@ref).
- `title`: defaults to the original's.

# Example

```jldoctest
julia> nom = additive_nomogram(
           Variable(:a; domain = (0, 10)),
           Variable(:b; domain = (0, 10)),
           Variable(:s),
       );

julia> skew = project(nom, [(0, 0), (70, 0), (70, 100), (0, 100)] =>
                           [(0, 0), (70, 20), (70, 80), (0, 100)]);

julia> round(evaluate(skew; a = 3, b = 4).s; digits = 10)
7.0
```
"""
function project(
    nom::Nomogram,
    H::AbstractMatrix;
    width::Real = nom.width,
    height::Real = nom.height,
    title::AbstractString = nom.title,
    preserve_aspect::Bool = false,
    nsamples::Int = 129,
)
    size(H) == (3, 3) ||
        throw(ArgumentError("a homography is a 3x3 matrix, got $(size(H))"))
    all(isfinite, H) || throw(ArgumentError("the homography has non-finite entries"))
    d = LinearAlgebra.det(H)
    abs(d) > 1e-12 || throw(
        ArgumentError(
            "the homography is singular (determinant $d), so it collapses the plane " *
            "rather than reshaping it",
        ),
    )
    width > 0 || throw(ArgumentError("width must be positive, got $width"))
    height > 0 || throw(ArgumentError("height must be positive, got $height"))

    vars = [s.variable for s in nom.scales]
    curves = Any[]
    for s in nom.scales
        n = name(s)
        f = s.curve
        push!(curves, t -> _project_point(H, _curvepoint(f, t), n))
    end
    scales, _ = _fitted_scales(
        vars,
        Tuple{Float64,Float64}[s.domain for s in nom.scales],
        curves,
        float(width),
        float(height),
        preserve_aspect,
        nsamples,
    )
    return Nomogram(
        scales;
        width = float(width),
        height = float(height),
        title,
        relation = nom.relation,
    )
end

function project(nom::Nomogram, quad::Pair; kwargs...)
    return project(nom, projective_transform(first(quad), last(quad)); kwargs...)
end
