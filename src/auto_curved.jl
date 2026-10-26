# --- charts with one curved scale -------------------------------------------------------
#
# A formula with no chart on three straight scales may still have one where two scales are
# straight and the third is a curve -- the roots of a quadratic, say. Write the relation as
# G(u, v, w) = 0, with u and v on the straight scales and w on the curve. Such a chart
# exists exactly when there are transforms f(u), g(v) under which every level set w = c is
# a straight line in the (f, g) plane, g = m(c) f + k(c); that is the classical Cauchy form,
# and the line for c meets the two straight scales -- f on x = 0, g on x = 1 -- in a
# family of lines all passing through (1/(1 - m), k/(1 - m)), which is where c goes.
#
# Along a level set, dg/df = -(g'(v)/f'(u)) G_u/G_v, so the lines are straight exactly when
# L = ln|G_u/G_v| splits as phi(w) + a(u) - b(v), with a = ln f', b = ln g', phi = ln|m|.
# Differentiating once in each of u and v gives, along each level set w = c,
#
#     L_uv = A w_u w_v + B w_uv,     A = phi''(c), B = phi'(c),
#
# so the test is a two-parameter least-squares fit on each of several level sets. The
# derivatives are third order in G and are taken exactly, by automatic differentiation of
# the formula as given. With B in hand, a'(u) = L_u - B w_u and b'(v) = -(L_v - B w_v),
# and two integrations give the straight scales.

const _FD = ForwardDiff

# The relation for one choice of curved variable, from z = F(x, y): the curved one last.
# Returns the generic G (for differentiation), the Float64 one, and which variables play u,
# v and w.
function _curved_roles(raw, F::_Formula, curved::Symbol)
    if curved === :z
        return ((u, v, w) -> w - raw(u, v), (u, v, w) -> w - F(u, v), (:x, :y, :z))
    elseif curved === :x
        return ((u, v, w) -> v - raw(w, u), (u, v, w) -> v - F(w, u), (:y, :z, :x))
    else
        return ((u, v, w) -> v - raw(u, w), (u, v, w) -> v - F(u, w), (:x, :z, :y))
    end
end

# w(u, v) near a point where w = c: Newton from c, in whatever number type (u, v) carry,
# so dual numbers come out with the derivatives of the implicit function. The branch is the
# one through c, which matters when w is double-valued -- the roots of a quadratic.
function _implicit_w(G, c::Float64)
    return x -> begin
        w = c + zero(x[1]) + zero(x[2])
        for _ = 1:8
            w -= G(x[1], x[2], w) / _FD.derivative(t -> G(x[1], x[2], t), w)
        end
        w
    end
end

function _logratio(G, W)
    return x -> begin
        w = W(x)
        gr = _FD.gradient(t -> G(t[1], t[2], t[3]), [x[1], x[2], w])
        log(abs(gr[1] / gr[2]))
    end
end

# A value of v on the level set w = c at u, or nothing.
function _level_v(Gf, u, c, vd)
    f(v) = Gf(u, v, c)
    a, b = f(vd[1]), f(vd[2])
    (isfinite(a) && isfinite(b)) || return nothing
    iszero(a) && return vd[1]
    iszero(b) && return vd[2]
    sign(a) == sign(b) && return nothing
    return Roots.find_zero(f, vd)
end

function _level_u(Gf, v, c, ud)
    roots = Roots.find_zeros(u -> Gf(u, v, c), ud...)
    return isempty(roots) ? nothing : roots[1]
end

_levels(lo, hi, n) =
    lo > 0 && hi / lo > 10 ? [lo * (hi / lo)^(i / (n - 1)) for i = 0:(n-1)] :
    collect(range(lo, hi; length = n))

# Tabulate a function known at points, with derivatives from central differences.
function _table(xs::Vector{Float64}, ys::Vector{Float64})
    n = length(xs)
    dys = similar(ys)
    for i = 1:n
        j, k = max(i - 1, 1), min(i + 1, n)
        dys[i] = (ys[k] - ys[j]) / (xs[k] - xs[j])
    end
    return TabulatedScale(xs, ys, dys)
end

# Build the chart for one choice of curved variable, or throw an ArgumentError saying in a
# line why it will not do.
function _curved_chart(
    F::_Formula,
    vx::Variable,
    vy::Variable,
    vz::Variable,
    zd::Tuple{Float64,Float64},
    curved::Symbol;
    width,
    height,
    title,
    rtol,
    reach,
    verify = true,
)
    G, Gf, (un, vn, wn) = _curved_roles(F.raw, F, curved)
    byname = Dict(:x => vx, :y => vy, :z => vz)
    doms = Dict(:x => vx.domain, :y => vy.domain, :z => zd)
    ud, vd, wd = doms[un], doms[vn], doms[wn]
    inset(d, f) = (d[1] + f * (d[2] - d[1]), d[2] - f * (d[2] - d[1]))

    # 1. The test: fit A, B on each level set.
    us = range(inset(ud, 0.01)...; length = 11)
    cs = _levels(inset(wd, 0.02)..., 25)
    fits = Tuple{Float64,Float64}[]
    for c in cs
        W = _implicit_w(G, c)
        L = _logratio(G, W)
        rows, rhs = Vector{Float64}[], Float64[]
        for u in us
            v = _level_v(Gf, u, c, vd)
            isnothing(v) && continue
            x = [u, v]
            push!(rhs, _FD.hessian(L, x)[1, 2])
            gw, hw = _FD.gradient(W, x), _FD.hessian(W, x)
            push!(rows, [gw[1] * gw[2], hw[1, 2]])
        end
        length(rhs) >= 4 || continue
        M = permutedims(reduce(hcat, rows))
        ab = M \ rhs
        scale = LinearAlgebra.norm(rhs)
        r = scale < 1e-12 ? 0.0 : LinearAlgebra.norm(M * ab - rhs) / scale
        isfinite(r) || throw(
            ArgumentError(
                "the formula is flat in :$un or :$vn somewhere, so they cannot both be " *
                "straight scales",
            ),
        )
        r <= 1e-6 || throw(
            ArgumentError(
                "its level sets are not straight lines under any transforms of " *
                ":$un and :$vn (misfit $(round(r; sigdigits = 2)))",
            ),
        )
        push!(fits, (c, ab[2]))
    end
    length(fits) >= 8 ||
        throw(ArgumentError("too few of its level sets cross the straight scales' domains"))
    B = _table(first.(fits), last.(fits))
    # Points for the scales are sought on levels across the whole range, not only those
    # fitted: near its ends a straight scale's values are reached only by levels the fit
    # left out, and B extends smoothly to them. Middle levels first.
    mid = (wd[1] + wd[2]) / 2
    levels = sort(_levels(wd..., 65); by = c -> abs(c - mid))

    # 2. The straight scales: log-slopes a' and b', integrated twice.
    # Values in the margin past what the formula reaches have no level set through them;
    # they take the slope at the nearest value that does.
    clampto(t, name) = clamp(t, reach[name]...)
    function aprime(u)
        u = clampto(u, un)
        for c in levels
            v = _level_v(Gf, u, c, vd)
            isnothing(v) && continue
            W = _implicit_w(G, c)
            x = [u, v]
            return _FD.gradient(_logratio(G, W), x)[1] -
                   transform(B, c) * _FD.gradient(W, x)[1]
        end
        throw(ArgumentError("no level set reaches :$un = $u"))
    end
    function bprime(v)
        v = clampto(v, vn)
        for c in levels
            u = _level_u(Gf, v, c, ud)
            isnothing(u) && continue
            W = _implicit_w(G, c)
            x = [u, v]
            return -(
                _FD.gradient(_logratio(G, W), x)[2] -
                transform(B, c) * _FD.gradient(W, x)[2]
            )
        end
        throw(ArgumentError("no level set reaches :$vn = $v"))
    end
    atab = _integrate_to_scale(aprime, ud...; n = 24)
    btab = _integrate_to_scale(bprime, vd...; n = 24)
    ftab = _integrate_to_scale(u -> exp(transform(atab, u)), ud...)
    gtab = _integrate_to_scale(v -> exp(transform(btab, v)), vd...)

    # 3. The curve: each level's line g = m f + k, from one point on it.
    mc, kc, ccs = Float64[], Float64[], Float64[]
    # A point on each level set. Near the ends of the curved variable's range its level set
    # shrinks into a corner of the other two domains, so look across them in full, ends
    # included, from either straight scale.
    ugrid = [(ud[1] + ud[2]) / 2; collect(range(ud...; length = 33))]
    vgrid = collect(range(vd...; length = 33))
    function onlevel(c)
        for u in ugrid
            v = _level_v(Gf, u, c, vd)
            isnothing(v) || return (u, v)
        end
        for v in vgrid
            u = _level_u(Gf, v, c, ud)
            isnothing(u) || return (u, v)
        end
        return nothing
    end
    # Levels for the curve run over exactly the values reached -- not the margin past them,
    # which no level set reaches -- with the ends themselves included, and closer together
    # toward them (Chebyshev spacing), where a level set shrinks into a corner and the
    # curve turns fastest. The table then interpolates everywhere a reading can land and
    # extrapolates only across the margin.
    rlo, rhi = reach[wn]
    for c in [(rlo + rhi) / 2 - (rhi - rlo) / 2 * cos(pi * k / 160) for k = 0:160]
        uv = onlevel(c)
        isnothing(uv) && continue
        u, v = uv
        gr = _FD.gradient(t -> G(t[1], t[2], t[3]), [u, v, c])
        m = -(exp(transform(btab, v)) / exp(transform(atab, u))) * gr[1] / gr[2]
        push!(mc, m)
        push!(kc, transform(gtab, v) - m * transform(ftab, u))
        push!(ccs, c)
    end
    length(ccs) >= 8 ||
        throw(ArgumentError("too few of its values could be placed on the curve"))
    # f and g are fixed only up to an affine map each, which moves the curve and not the
    # answers. Map each straight scale onto [0, 1], so both run the full height, and orient
    # g so the curve lies between them where the lines allow: m < 0 does that. A line
    # g = m f + k becomes g~ = (s m / df) f~ + s (k + m f0 - g0) / dg, with f~ = (f - f0)/df
    # and g~ = s (g - g0)/dg.
    affine(t, a, b) = TabulatedScale(t.x, a .* t.y .+ b, a .* t.dy)
    f0, df = transform(ftab, ud[1]), transform(ftab, ud[2]) - transform(ftab, ud[1])
    g0, dg = transform(gtab, vd[1]), transform(gtab, vd[2]) - transform(gtab, vd[1])
    s = count(<(0), mc .* df ./ dg) >= length(mc) / 2 ? 1.0 : -1.0
    if s < 0
        g0, dg = g0 + dg, -dg
    end
    kc = (kc .+ mc .* f0 .- g0) ./ dg
    mc = mc .* df ./ dg
    ftab = affine(ftab, 1 / df, -f0 / df)
    gtab = affine(gtab, 1 / dg, -g0 / dg)
    d = 1 .- mc
    (all(>(0), d) || all(<(0), d)) && minimum(abs, d) > 1e-6 || throw(
        ArgumentError("its curve would run off to infinity between two of its values"),
    )
    mtab, ktab = _table(ccs, mc), _table(ccs, kc)
    curve(c) = (m = transform(mtab, c); (1 / (1 - m), transform(ktab, c) / (1 - m)))

    # 4. The chart, and the check every automatic chart passes.
    vw = byname[wn]
    wstyle =
        vw.tickstyle !== :auto ? vw.tickstyle :
        wd[1] > 0 && wd[2] / wd[1] >= 10 ? :log : :linear
    nom = determinant_nomogram(
        _with_derived(byname[un], ftab, ud) => u -> (0.0, transform(ftab, u)),
        _with_derived(byname[vn], gtab, vd) => v -> (1.0, transform(gtab, v)),
        rebuild(vw; domain = wd, tickstyle = wstyle) => curve;
        width,
        height,
        title,
        relation = Gf,
    )
    verify || return nom
    zspan = zd[2] - zd[1]
    worst = 0.0
    for a in range(vx.domain...; length = 9), b in range(vy.domain...; length = 9)
        got = try
            _evaluate(nom; NamedTuple{(vx.name, vy.name)}((a, b))...)[vz.name]
        catch err
            err isa ArgumentError || rethrow()
            throw(
                ArgumentError(
                    "the chart built from it cannot be read at $(vx.name) = $a, " *
                    "$(vy.name) = $b: $(first(split(err.msg, '.')))",
                ),
            )
        end
        worst = max(worst, abs(got - F(a, b)) / zspan)
    end
    worst <= rtol || throw(
        ArgumentError(
            "the chart built from it disagrees with F by $(round(worst; sigdigits = 3)) " *
            "of the result's range",
        ),
    )
    return nom
end

# The range of F over the input rectangle, which the result's scale has to cover. With
# a curved input F need not be monotone, so an extreme can lie inside the rectangle: each
# is narrowed down by sampling ever smaller boxes around the best point found.
function _value_range(F, xd, yd; n = 41)
    function extreme(sgn)
        best(xs, ys) = argmax(p -> sgn * F(p...), ((a, b) for a in xs, b in ys))
        p = best(range(xd...; length = n), range(yd...; length = n))
        hx, hy = (xd[2] - xd[1]) / (n - 1), (yd[2] - yd[1]) / (n - 1)
        for _ = 1:6
            p = best(
                range(max(xd[1], p[1] - hx), min(xd[2], p[1] + hx); length = 11),
                range(max(yd[1], p[2] - hy), min(yd[2], p[2] + hy); length = 11),
            )
            hx, hy = hx / 5, hy / 5
        end
        return F(p...)
    end
    return (extreme(-1), extreme(1))
end

# Every choice of curved variable, result first; the first that gives a verified chart.
function _auto_curved(F::_Formula, vx, vy, vz; width, height, title, rtol, straight)
    zd = _value_range(F, vx.domain, vy.domain)
    if !isnothing(vz.domain)
        span = zd[2] - zd[1]
        ok = all(i -> abs(vz.domain[i] - zd[i]) <= 1e-6 * span, 1:2)
        ok || throw(
            ArgumentError(
                "domain $(vz.domain) for :$(vz.name) does not match what F reaches over " *
                "the input domains, ($(default_tickformat(zd[1])), " *
                "$(default_tickformat(zd[2]))). Leave it unset to have it derived.",
            ),
        )
        zd = vz.domain
    end
    reach = Dict(:x => vx.domain, :y => vy.domain, :z => zd)
    # A margin past the values the formula reaches: the scales are tabulated, accurate to
    # around 1e-7, and a reading at a corner of the inputs -- exactly the end of the range
    # -- must not land a hair outside the result's scale and be refused.
    pad = 1e-5 * (zd[2] - zd[1])
    zd = (zd[1] - pad, zd[2] + pad)
    reasons = String[]
    for curved in (:z, :x, :y)
        cname = curved === :z ? vz.name : curved === :x ? vx.name : vy.name
        try
            return _curved_chart(
                F,
                vx,
                vy,
                vz,
                zd,
                curved;
                width,
                height,
                title,
                rtol,
                reach,
            )
        catch err
            msg = if err isa ArgumentError
                err.msg
            elseif (err isa MethodError || err isa InexactError) &&
                   occursin("Dual", sprint(showerror, err))
                "the formula could not be differentiated automatically; it has to " *
                "accept any number type, not only Float64"
            elseif err isa DomainError
                "the formula is undefined somewhere it had to be evaluated"
            else
                rethrow()
            end
            push!(reasons, "  • :$cname on a curve — $msg")
        end
    end
    throw(
        ArgumentError(
            "$(vz.name) = F($(vx.name), $(vy.name)) has no chart on three straight " *
            "scales ($straight), nor one with a single curved scale:\n" *
            join(reasons, "\n") *
            "\nIt may still have a chart with more curved scales, or a grid, from a form " *
            "you derive: see `determinant_nomogram` and `grid_nomogram`. Or draw it as " *
            "curves of constant value, which any formula has: see `contour_chart`.",
        ),
    )
end
