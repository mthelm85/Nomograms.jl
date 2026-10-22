"""
    ScaleTransform

A strictly monotone map from a variable's own units onto the *functional* axis of a
nomogram, together with its inverse.

A three-variable alignment chart exists whenever the relation can be written in the
separable form `f1(u) + f2(v) = f3(w)`; a `ScaleTransform` is one of those `f`s. Bundling
the inverse with the forward map lets the construction derive the third variable's domain
exactly instead of searching for it.

Implement [`transform`](@ref) for a new subtype, and either [`invtransform`](@ref) or
`hasinverse(::T) = false` (in which case inversion falls back to bracketed root finding
over the variable's domain).
"""
abstract type ScaleTransform end

"""
    transform(t::ScaleTransform, x)

Map a value in the variable's own units onto the functional axis.
"""
function transform end

"""
    invtransform(t::ScaleTransform, y)

Map a functional-axis value back to the variable's own units. Only defined for transforms
with `hasinverse(t) == true`.
"""
function invtransform end

"""
    hasinverse(t::ScaleTransform) -> Bool

Whether `t` provides a closed-form [`invtransform`](@ref).
"""
hasinverse(::ScaleTransform) = true

"""
    tickstyle(t::ScaleTransform) -> Symbol

How tick values should be chosen for a variable carrying this transform: `:linear` picks
evenly spaced "nice" numbers, `:log` picks decades subdivided by mantissa.
"""
tickstyle(::ScaleTransform) = :linear

"""
    LinearScale(a = 1, b = 0)

`f(x) = a*x + b`. The default transform, used for plain additive relations.
"""
struct LinearScale <: ScaleTransform
    a::Float64
    b::Float64
end
LinearScale(; a = 1.0, b = 0.0) = LinearScale(a, b)
LinearScale(a) = LinearScale(a, 0.0)
transform(t::LinearScale, x) = t.a * x + t.b
invtransform(t::LinearScale, y) = (y - t.b) / t.a

"""
    LogScale(base = 10)

`f(x) = log_base(x)`, defined for `x > 0`. Turns a product relation into an additive one.
"""
struct LogScale <: ScaleTransform
    base::Float64
end
LogScale() = LogScale(10.0)
transform(t::LogScale, x) = log(x) / log(t.base)
invtransform(t::LogScale, y) = t.base^y
tickstyle(::LogScale) = :log

"""
    ReciprocalScale()

`f(x) = 1/x`, defined away from zero. The transform behind parallel-resistance and
thin-lens style charts.
"""
struct ReciprocalScale <: ScaleTransform end
transform(::ReciprocalScale, x) = inv(x)
invtransform(::ReciprocalScale, y) = inv(y)

"""
    PowerScale(p)

`f(x) = x^p`.
"""
struct PowerScale <: ScaleTransform
    p::Float64
end
transform(t::PowerScale, x) = x^t.p
invtransform(t::PowerScale, y) = y^inv(t.p)

"""
    FunctionScale(f, finv = nothing)

An arbitrary monotone `f` supplied by the user. Without `finv`, inversion is done by
bracketed root finding over the variable's domain, so such a variable must declare its
`domain` explicitly.
"""
struct FunctionScale{F,G} <: ScaleTransform
    f::F
    finv::G
end
FunctionScale(f) = FunctionScale(f, nothing)
transform(t::FunctionScale, x) = t.f(x)
invtransform(t::FunctionScale, y) = t.finv(y)
hasinverse(::FunctionScale{F,Nothing}) where {F} = false

"""
    ComposedScale(outer, inner)

`f(x) = outer(inner(x))`. Used internally by [`multiplicative_nomogram`](@ref) to wrap a
user transform in a logarithm, and available for building transforms such as
`ComposedScale(LogScale(), PowerScale(2))`.
"""
struct ComposedScale{O<:ScaleTransform,I<:ScaleTransform} <: ScaleTransform
    outer::O
    inner::I
end
transform(t::ComposedScale, x) = transform(t.outer, transform(t.inner, x))
invtransform(t::ComposedScale, y) = invtransform(t.inner, invtransform(t.outer, y))
hasinverse(t::ComposedScale) = hasinverse(t.outer) && hasinverse(t.inner)
tickstyle(t::ComposedScale) = tickstyle(t.outer)

"""
    invert(t::ScaleTransform, y, domain) -> Real

Solve `transform(t, x) == y` for `x`. Uses the closed-form inverse when available and
falls back to bracketed root finding over `domain = (lo, hi)` otherwise.
"""
function invert(t::ScaleTransform, y, domain::Tuple{Real,Real})
    hasinverse(t) && return float(invtransform(t, y))
    lo, hi = float(domain[1]), float(domain[2])
    g = x -> transform(t, x) - y
    glo, ghi = g(lo), g(hi)
    if glo * ghi > 0
        throw(
            ArgumentError(
                "cannot invert transform at $y: no sign change over domain " *
                "($lo, $hi) (residuals $glo and $ghi). Supply an inverse via " *
                "FunctionScale(f, finv) or widen the domain.",
            ),
        )
    end
    return Roots.find_zero(g, (lo, hi), Roots.A42())
end

"""
    frange(t::ScaleTransform, domain) -> (lo, hi)

The interval spanned by `transform(t, .)` over `domain`, after checking that `t` really is
monotone there. A fold-over would make the axis double back on itself, which no amount of
drawing can rescue, so it is rejected early.
"""
function frange(t::ScaleTransform, domain::Tuple{Real,Real}; nsamples::Int = 65)
    lo, hi = float(domain[1]), float(domain[2])
    lo < hi || throw(ArgumentError("domain must be increasing, got ($lo, $hi)"))
    xs = range(lo, hi; length = nsamples)
    ys = [float(transform(t, x)) for x in xs]
    all(isfinite, ys) || throw(
        ArgumentError(
            "transform $(typeof(t)) is not finite everywhere on ($lo, $hi); " *
            "check that the domain excludes singularities",
        ),
    )
    d = diff(ys)
    if !(all(>(0), d) || all(<(0), d))
        throw(
            ArgumentError(
                "transform $(typeof(t)) is not monotone on ($lo, $hi); a nomogram " *
                "axis must be a one-to-one map from values to positions",
            ),
        )
    end
    return (min(ys[1], ys[end]), max(ys[1], ys[end]))
end
