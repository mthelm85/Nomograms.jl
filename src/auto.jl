"""
    TabulatedScale

A [`ScaleTransform`](@ref) known only numerically: values and derivatives at a set of
nodes, joined by cubic Hermite interpolation.

This is what [`auto_nomogram`](@ref) builds when it derives a transform by integration,
since most of those integrals have no closed form worth looking for. With exact derivatives
at the nodes the interpolation error is fourth order in the node spacing, far below
anything a chart can show. Outside the tabulated range it extrapolates linearly, so a
solver probing an ulp past a domain edge gets an answer rather than an exception.
"""
struct TabulatedScale <: ScaleTransform
    x::Vector{Float64}
    y::Vector{Float64}
    dy::Vector{Float64}
end

function _hermite(x0, x1, y0, y1, d0, d1, x)
    h = x1 - x0
    t = (x - x0) / h
    t2 = t * t
    t3 = t2 * t
    return (2t3 - 3t2 + 1) * y0 +
           (t3 - 2t2 + t) * h * d0 +
           (-2t3 + 3t2) * y1 +
           (t3 - t2) * h * d1
end

function transform(t::TabulatedScale, x)
    x = float(x)
    xs = t.x
    x <= xs[1] && return t.y[1] + t.dy[1] * (x - xs[1])
    x >= xs[end] && return t.y[end] + t.dy[end] * (x - xs[end])
    i = clamp(searchsortedlast(xs, x), 1, length(xs) - 1)
    return _hermite(xs[i], xs[i+1], t.y[i], t.y[i+1], t.dy[i], t.dy[i+1], x)
end

function invtransform(t::TabulatedScale, y)
    y = float(y)
    ys = t.y
    up = ys[end] > ys[1]
    if up ? y <= ys[1] : y >= ys[1]
        return t.x[1] + (y - ys[1]) / t.dy[1]
    elseif up ? y >= ys[end] : y <= ys[end]
        return t.x[end] + (y - ys[end]) / t.dy[end]
    end
    i = up ? searchsortedlast(ys, y) : searchsortedlast(ys, y; rev = true)
    n = length(ys)
    # A neighbouring cell either side, since a cubic can overshoot its end values by a
    # whisker and leave the root just across a node.
    lo, hi = t.x[max(i - 1, 1)], t.x[min(i + 2, n)]
    return Roots.find_zero(v -> transform(t, v) - y, (lo, hi), Roots.A42())
end

function Base.show(io::IO, t::TabulatedScale)
    print(io, "TabulatedScale(", length(t.x), " nodes on [")
    print(io, default_tickformat(t.x[1]), ", ", default_tickformat(t.x[end]), "])")
    return nothing
end

# --- numerical calculus ---------------------------------------------------------------

# Derivative by central differences, with a second-order one-sided formula at the edges,
# so that a function is never evaluated outside the region its caller vouched for. A
# formula like log(r) may simply not exist past the end of its domain.
function _deriv(f, x::Float64, lo::Float64, hi::Float64)
    # The step scales with x itself, not only the domain. A span-relative step is a large
    # fraction of x at the small end of a domain spanning decades -- 1% of R at R = 1 over
    # [1, 1000] -- and the truncation error there was enough for V²/R to disagree with its
    # own chart. The floor keeps the step sensible where x passes through zero.
    h = 6e-6 * max(abs(x), 1e-2 * (hi - lo))
    if x - h >= lo && x + h <= hi
        return (f(x + h) - f(x - h)) / (2h)
    elseif x + 2h <= hi
        return (-3f(x) + 4f(x + h) - f(x + 2h)) / (2h)
    else
        return (3f(x) - 4f(x - h) + f(x - 2h)) / (2h)
    end
end

const _GL5 = (
    (-0.9061798459386640, 0.2369268850561891),
    (-0.5384693101056831, 0.4786286704993665),
    (0.0, 0.5688888888888889),
    (0.5384693101056831, 0.4786286704993665),
    (0.9061798459386640, 0.2369268850561891),
)

function _gauss5(f, a::Float64, b::Float64)
    m, r = (a + b) / 2, (b - a) / 2
    return r * sum(w * f(m + r * t) for (t, w) in _GL5)
end

# Nodes for tabulating a transform. Geometric when the domain spans a decade or more, since
# that is where a log-like integrand changes fastest and uniform panels would waste most of
# their points on the flat end.
function _nodes(lo::Float64, hi::Float64, n::Int)
    xs = if lo > 0 && hi / lo > 10
        [lo * (hi / lo)^(i / n) for i = 0:n]
    else
        collect(range(lo, hi; length = n + 1))
    end
    xs[1], xs[end] = lo, hi
    return xs
end

# Integrate a derivative into a transform: f(lo) = 0, and f(x) = ∫ f′ from lo.
function _integrate_to_scale(deriv, lo::Float64, hi::Float64; n::Int = 256)
    xs = _nodes(lo, hi, n)
    ys = zeros(n + 1)
    for i = 1:n
        ys[i+1] = ys[i] + _gauss5(deriv, xs[i], xs[i+1])
    end
    return TabulatedScale(xs, ys, [deriv(x) for x in xs])
end

# --- the criterion --------------------------------------------------------------------

struct _Partials{F}
    F::F
    xd::Tuple{Float64,Float64}
    yd::Tuple{Float64,Float64}
end

_Fx(P::_Partials, x, y) = _deriv(t -> P.F(t, y), float(x), P.xd...)
_Fy(P::_Partials, x, y) = _deriv(t -> P.F(x, t), float(y), P.yd...)

function _not_monotone(name::Symbol, dom)
    return ArgumentError(
        "the function is not monotone in :$name over $dom, so :$name cannot be a " *
        "straight scale: somewhere along it, increasing :$name stops moving the result " *
        "the same way. Narrow the domain to a range where it is monotone.",
    )
end

"""
    separability_defect(F, xdomain, ydomain; n = 9) -> Float64

How far `z = F(x, y)` is from being drawable on three straight parallel scales. Zero (to
rounding) when transforms `f`, `g`, `h` exist with `h(z) = f(x) + g(y)`; otherwise the
largest violation found on an `n × n` grid.

The test is the classical one: such transforms exist exactly when `ln|F_x / F_y|` splits
into a function of `x` plus a function of `y`. Rather than take a mixed second derivative
of that — numerically noisy — this checks the equivalent identity

```
L(x, y) − L(x, y₀) − L(x₀, y) + L(x₀, y₀) = 0
```

which needs only first derivatives. Throws if `F` is not monotone in each variable, since
then no straight scale can carry it whatever the transforms.
"""
function separability_defect(F, xd, yd; n::Int = 9, xname = :x, yname = :y)
    xd = (float(xd[1]), float(xd[2]))
    yd = (float(yd[1]), float(yd[2]))
    P = _Partials((a, b) -> float(F(a, b)), xd, yd)
    inset(d) = range(d[1] + 0.02 * (d[2] - d[1]), d[2] - 0.02 * (d[2] - d[1]); length = n)
    xs, ys = inset(xd), inset(yd)

    fx = [_Fx(P, x, y) for x in xs, y in ys]
    fy = [_Fy(P, x, y) for x in xs, y in ys]
    all(isfinite, fx) && all(isfinite, fy) || throw(
        ArgumentError(
            "the function is not finite everywhere on the domains given; check them " *
            "for singularities",
        ),
    )
    for (d, nm, dom) in ((fx, xname, xd), (fy, yname, yd))
        if all(iszero, d)
            throw(
                ArgumentError(
                    "the result does not depend on :$nm at all, so there is nothing " *
                    "for its scale to show",
                ),
            )
        end
        (all(>(0), d) || all(<(0), d)) || throw(_not_monotone(nm, dom))
    end

    L = log.(abs.(fx ./ fy))
    r = (n + 1) ÷ 2
    worst = 0.0
    for i = 1:n, j = 1:n
        worst = max(worst, abs(L[i, j] - L[i, r] - L[r, j] + L[r, r]))
    end
    return worst
end

# --- building it ----------------------------------------------------------------------

# The user's formula behind one concrete type per arity. Every closure is a type of its
# own, so code specialised on the formula would be compiled afresh for each one a user
# writes, which is most of the wait on a first call. Behind this wrapper it compiles once
# and is cached at precompile time. A FunctionWrapper calls through a pointer of a fixed
# signature, so unlike a field typed `Any` it costs no dynamic dispatch and no boxed
# result -- and the formula is called a great many times. The stages a chain builds
# internally stay plain closures: their types depend on the shape of the chain, not on the
# formula, so the cache covers them too.
struct _Formula{N}
    f::FunctionWrapper{Float64,NTuple{N,Float64}}
end
_Formula(f::_Formula, n::Int) = f
_Formula(f, n::Int) = _Formula{n}(
    FunctionWrapper{Float64,NTuple{n,Float64}}((args...) -> Float64(f(args...))),
)
(g::_Formula{N})(args::Vararg{Float64,N}) where {N} = g.f(args...)
(g::_Formula{N})(args::Vararg{Real,N}) where {N} = g.f(map(Float64, args)...)

_as_variable(v::Variable) = v
_as_variable(p::Pair{Symbol}) = Variable(first(p); domain = last(p))
_as_variable(s::Symbol) = Variable(s)

# Which kind of round number suits an axis. Two questions, in order.
#
# Is it shaped like a log axis? If the transform puts the geometric mean of the domain
# nearer the middle of the axis than the arithmetic mean, it spreads values multiplicatively.
# This used to be asked only of domains spanning a decade, which left log(b, x)'s base axis
# over [2, 10] -- strongly non-linear, being ln ln b -- on linear ticks with nothing
# labelled between 2 and 4.
#
# And if so, does it actually label better? Shape alone is not enough: over 100 to 200 the
# whole log mantissas are only the two ends. So both styles are thinned onto a nominal axis
# and the one that lands more labels wins, ties going to the log.
# Whether a transform spreads a domain logarithmically: the geometric mean lands nearer
# the middle of the axis than the arithmetic mean does.
function _log_shaped(t::ScaleTransform, dom::Tuple{Float64,Float64})
    lo, hi = dom
    lo > 0 || return false
    f0, f1 = transform(t, lo), transform(t, hi)
    frac(v) = (transform(t, v) - f0) / (f1 - f0)
    return abs(frac(sqrt(lo * hi)) - 0.5) < abs(frac((lo + hi) / 2) - 0.5)
end

function _auto_tickstyle(t::ScaleTransform, dom::Tuple{Float64,Float64})
    lo, hi = dom
    _log_shaped(t, dom) || return :linear
    f0, f1 = transform(t, lo), transform(t, hi)
    frac(v) = (transform(t, v) - f0) / (f1 - f0)
    probe = Scale(Variable(:_probe; domain = dom), dom, v -> (0.0, 100 * frac(v)), :left)
    labelled((majors, minors)) =
        count(k -> !isnothing(k.label), thin_ticks(probe, majors, minors, NomogramStyle()))
    return labelled(log_ticks(lo, hi)) >= labelled(linear_ticks(lo, hi)) ? :log : :linear
end

function _with_derived(v::Variable, t::ScaleTransform, dom)
    style =
        (v.tickstyle === :auto && isnothing(v.ticks)) ? _auto_tickstyle(t, dom) :
        v.tickstyle
    return rebuild(v; transform = t, domain = dom, tickstyle = style)
end

# A straight path across the domain rectangle along which F rises steadily: from the
# corner minimising F to the corner maximising it, so every step moves both inputs in the
# direction that raises F, and F is monotone in the path parameter whatever its shape. The
# parameter is solved for at nodes placed in F's own values -- geometrically when F spans
# decades, since nodes spaced evenly along the path leave the small end of a wide-ranging
# result (x*y over [1, 100]²) with too few to hold the accuracy.
function _rising_path(Fc, P::_Partials, xd, yd)
    xr, yr = (xd[1] + xd[2]) / 2, (yd[1] + yd[2]) / 2
    sx = sign(_Fx(P, xr, yr))
    sy = sign(_Fy(P, xr, yr))
    cmin = (sx > 0 ? xd[1] : xd[2], sy > 0 ? yd[1] : yd[2])
    cmax = (sx > 0 ? xd[2] : xd[1], sy > 0 ? yd[2] : yd[1])
    dx, dy = cmax[1] - cmin[1], cmax[2] - cmin[2]
    along(s) = (cmin[1] + s * dx, cmin[2] + s * dy)
    zlo, zhi = Fc(cmin...), Fc(cmax...)
    zlo < zhi || throw(
        ArgumentError(
            "the result does not increase across the domain, which a separable, " *
            "monotone function must; the domains may straddle a singularity",
        ),
    )
    zs = _nodes(zlo, zhi, 1024)
    svals = [
        i == 1 ? 0.0 :
        i == length(zs) ? 1.0 :
        Roots.find_zero(s -> Fc(along(s)...) - zs[i], (0.0, 1.0), Roots.A42()) for
        i in eachindex(zs)
    ]
    return (; along, zs, svals, dx, dy)
end

"""
    auto_nomogram(F, x, y, [z]; kwargs...) -> Nomogram
    auto_nomogram(F, x, y, w, ..., z; kwargs...) -> CompoundNomogram

A nomogram for `z = F(x, y)` — or `z = F(x, y, w, …)` with any number of inputs — from the
function as you would ordinarily write it.

Every other constructor asks you to rewrite your formula first — into `f₁(u) + f₂(v) =
f₃(w)`, with the transforms chosen and the inverses supplied. This one does that part:

```julia
nom = auto_nomogram((r, n) -> (1 + r/100)^n, :r => (1, 20), :n => (1, 30), :growth)
```

Each input is a [`Variable`](@ref) with a domain, or simply `name => (lo, hi)`, in the order
`F` takes them. The last argument names the result — a `Variable` or a name — except that
with exactly two arguments after `F` both are inputs and the result is `:z`. The result's
domain is derived from `F` and, if you give one, checked. Labels, tick formats and explicit
ticks on the variables are kept. Transforms are not — deriving them is the point.

# How

`F` can be drawn on three straight parallel scales exactly when transforms exist with
`h(z) = f(x) + g(y)`. Differentiating shows that requires

```
F_x / F_y = f′(x) / g′(y)
```

so the test is whether `ln|F_x / F_y|` splits into a part in `x` plus a part in `y` (see
[`separability_defect`](@ref)). When it does, that same identity hands back `f′` and `g′`;
integrating them gives `f` and `g`, and `h` follows by walking the level sets of `F`. The
transforms come out as [`TabulatedScale`](@ref)s, and each axis gets linear or
decade-based ticks according to how its transform spreads its values.

# Verification

Whatever the numerics did, the finished chart is then *solved* over a grid of inputs and
every answer compared with `F` itself. A chart that disagrees by more than `rtol` of the
result's range is never returned — the call throws instead, quoting the error it measured.
The chart's `relation` is `F`, so [`relation_residual`](@ref) and
[`max_relation_residual`](@ref) check against your formula rather than the derived one.

# More than two inputs

One alignment line relates three scales, so `F(x, y, w)` is drawn as a chain of
[`compound`](@ref) stages and the result is a [`CompoundNomogram`](@ref). That needs `F` to
depend on some pair of its inputs only through a single combination of them, `t = T(x, y)`,
which holds exactly when `F_x / F_y` does not depend on `w`. The first stage then charts
`t`, the second charts `z` from `t` and `w`, and the reader pivots between them on a
bare turning line. More inputs are reduced the same way, a pair at a time.

Every pairing is tried. One can fail for two reasons: `F` may not depend on the pair
through one combination at all, or it may, but with the two stages wanting the
intermediate on different scales — `(x + y)·w` is the standard case, linear in the first
stage and logarithmic in the second. That second kind needs a net, which has to be derived
by hand ([`net_nomogram`](@ref)). If no pairing works, the error lists what went wrong
with each.

# When it can't

A formula that is not separable has no three-straight-scale chart at all, and the error
says so. It may still have one with a curved scale or a net, but those need a form derived
by hand: see [`determinant_nomogram`](@ref) and [`net_nomogram`](@ref). A formula that is
not monotone in a variable over its domain fails too, since a straight scale cannot fold
back on itself.

# Keywords
- `width`, `height`, `title`: as for the other constructions. `width` defaults to 70 for
  two inputs. In a chain it is the width of the last stage, and the default grows with
  the chain, 35 per input after the first. An earlier stage always has one input inside
  the rest of the chain, so it is mirrored and resized, between a quarter of `width` and
  `width`, to land that input mid-way between the axes already there. Every way of
  pairing the inputs is tried for the first two levels of the chain, and the one whose
  axes end up furthest apart is kept.
- `rtol = 1e-6`: the largest disagreement with `F` to accept, relative to the result's
  range.
"""
function auto_nomogram(
    F,
    args...;
    width::Union{Nothing,Real} = nothing,
    height::Real = 100.0,
    title::AbstractString = "",
    rtol::Real = 1e-6,
)
    length(args) >= 2 ||
        throw(ArgumentError("auto_nomogram needs at least two inputs, got $(length(args))"))
    inputs, result = length(args) == 2 ? (args, :z) : (args[1:(end-1)], args[end])
    vs = Variable[_as_variable(a) for a in inputs]
    vz = _as_variable(result)
    for v in vs
        isnothing(v.domain) &&
            throw(ArgumentError("input variable :$(v.name) must declare a domain"))
    end
    F, vz = _with_dimensions(F, vs, vz)
    F = _Formula(F, length(vs))
    names = [v.name for v in vs]
    allunique([names; vz.name]) ||
        throw(ArgumentError("variable names must be distinct, got $([names; vz.name])"))
    width = float(something(width, 35.0 * max(length(vs) - 1, 2)))
    length(vs) == 2 && return _auto_two(F, vs[1], vs[2], vz; width, height, title, rtol)
    return _auto_chain(F, vs, vz; width, height, title, rtol)
end

# With units, the formula is called with quantities -- so it may compute in whatever
# units it likes -- and its result converted to the result variable's unit, which is
# taken from the formula's own result if not given. That conversion is a dimension check
# on the formula, made once here, before anything is built; everything after works on
# plain numbers in each variable's unit.
function _with_dimensions(F, vs::Vector{Variable}, vz::Variable)
    all(v -> isnothing(v.unit), vs) && isnothing(vz.unit) && return (F, vz)
    Fq = (xs...) -> F((_attach(v.unit, x) for (v, x) in zip(vs, xs))...)
    out = Fq(((v.domain[1] + v.domain[2]) / 2 for v in vs)...)
    isnothing(vz.unit) && (vz = rebuild(vz; unit = _infer_unit(out)))
    try
        _strip(vz.unit, out, vz.name)
    catch err
        err isa ArgumentError || rethrow()
        names = join(string.(getfield.(vs, :name)), ", ")
        throw(
            ArgumentError(
                "$(vz.name) = F($names) does not come out in the units of " *
                ":$(vz.name): $(err.msg)",
            ),
        )
    end
    return ((xs...) -> _strip(vz.unit, Fq(xs...), vz.name), vz)
end

function _auto_two(
    F,
    vx::Variable,
    vy::Variable,
    vz::Variable;
    width::Real = 70.0,
    height::Real = 100.0,
    title::AbstractString = "",
    rtol::Real = 1e-6,
    verify::Bool = true,
)
    xd, yd = vx.domain, vy.domain
    Fc = (a, b) -> float(F(a, b))

    defect = separability_defect(Fc, xd, yd; xname = vx.name, yname = vy.name)
    if defect > 1e-6
        throw(
            ArgumentError(
                "$(vz.name) = F($(vx.name), $(vy.name)) is not separable: there are no " *
                "transforms with h($(vz.name)) = f($(vx.name)) + g($(vy.name)), so it " *
                "has no chart on three straight scales (ln|F_x/F_y| fails to split by " *
                "$(round(defect; sigdigits = 3))). It may still have one with a curved " *
                "scale or a net, which need a form you derive: see " *
                "`determinant_nomogram` and `net_nomogram`.",
            ),
        )
    end

    P = _Partials(Fc, xd, yd)
    xr, yr = (xd[1] + xd[2]) / 2, (yd[1] + yd[2]) / 2
    R(a, b) = _Fx(P, a, b) / _Fy(P, a, b)
    Rref = R(xr, yr)
    # F_x/F_y = f′(x)/g′(y), and it separates, so fixing one argument at a time gives
    # each derivative up to a common factor -- which is all a nomogram needs.
    fprime(a) = R(a, yr)
    gprime(b) = Rref / R(xr, b)
    tf = _integrate_to_scale(fprime, xd...)
    tg = _integrate_to_scale(gprime, yd...)

    # h along the rising path. Moving along it changes f and g the same way (h′F_x = f′
    # and h′F_y = g′ fix their signs together), so h is tabulated without ever having to
    # invert F.
    path = _rising_path(Fc, P, xd, yd)
    zs, svals, dx, dy = path.zs, path.svals, path.dx, path.dy
    zlo, zhi = zs[1], zs[end]
    pts = path.along.(svals)
    hs = [transform(tf, p[1]) + transform(tg, p[2]) for p in pts]
    hd = [
        (fprime(p[1]) * dx + gprime(p[2]) * dy) / (_Fx(P, p...) * dx + _Fy(P, p...) * dy) for p in pts
    ]
    th = TabulatedScale(zs, hs, hd)
    zd = (zlo, zhi)

    # The transforms are only fixed up to a common sign -- -f, -g, -h serve as well as f,
    # g, h -- and as derived, h decreases whenever F decreases in its second argument,
    # turning the result axis upside down: profit = revenue - cost, atan(y, x). Choose the
    # sign that makes the result read upwards.
    if hs[end] < hs[1]
        flip(t::TabulatedScale) = TabulatedScale(t.x, -t.y, -t.dy)
        tf, tg, th = flip(tf), flip(tg), flip(th)
    end

    if !isnothing(vz.domain)
        ok =
            isapprox(vz.domain[1], zd[1]; rtol = 1e-6, atol = 1e-12) &&
            isapprox(vz.domain[2], zd[2]; rtol = 1e-6, atol = 1e-12)
        ok || throw(
            ArgumentError(
                "domain $(vz.domain) for :$(vz.name) does not match what F reaches over " *
                "the input domains, ($(default_tickformat(zd[1])), " *
                "$(default_tickformat(zd[2]))). Leave it unset to have it derived.",
            ),
        )
    end

    built = additive_nomogram(
        _with_derived(vx, tf, xd),
        _with_derived(vy, tg, yd),
        _with_derived(vz, th, zd);
        width,
        height,
        title,
    )
    nom = Nomogram(
        built.scales;
        width = built.width,
        height = built.height,
        title = built.title,
        relation = (a, b, c) -> c - Fc(a, b),
    )

    # The chart is only worth returning if reading it gives F. Solve it and check -- unless
    # it is one stage of a chain, which is checked whole once assembled.
    verify || return nom
    zspan = zd[2] - zd[1]
    worst = 0.0
    names = (vx.name, vy.name)
    for a in range(xd...; length = 9), b in range(yd...; length = 9)
        got = _evaluate(nom; NamedTuple{names}((a, b))...)[vz.name]
        worst = max(worst, abs(got - Fc(a, b)) / zspan)
    end
    worst <= rtol || throw(
        ArgumentError(
            "the derived chart disagrees with F by $(round(worst; sigdigits = 3)) of the " *
            "result's range, more than rtol = $rtol, so it has not been returned. This " *
            "should not happen for a separable, smooth function; please report it.",
        ),
    )
    return nom
end
