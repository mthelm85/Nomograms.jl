# auto_nomogram for more than two inputs: split the formula into a chain of two-input
# charts linked by pivots, and hand each link to the two-input machinery.

function _partial_n(F, vals::Vector{Float64}, k::Int, dom::Tuple{Float64,Float64})
    return _deriv(t -> (v = copy(vals); v[k] = t; float(F(v...))), vals[k], dom...)
end

# How far F is from depending on inputs i and j only through one combination of them,
# F = G(T(xᵢ, xⱼ), others…). That holds exactly when the ratio F_i/F_j is the same whatever
# the other inputs are, because then F's level curves in the (xᵢ, xⱼ) plane are the same
# family of curves for every setting of the rest -- the level curves of T. Measured, like
# separability_defect, on the log of the ratio, so it is scale-free.
function _pivot_defect(
    F,
    doms::Vector{Tuple{Float64,Float64}},
    names::Vector{Symbol},
    i::Int,
    j::Int;
    n::Int = 5,
)
    m = length(doms)
    inset(d, k) =
        range(d[1] + 0.02 * (d[2] - d[1]), d[2] - 0.02 * (d[2] - d[1]); length = k)
    others = [k for k = 1:m if k != i && k != j]
    ref = [(d[1] + d[2]) / 2 for d in doms]
    restvalues = collect(Iterators.product((inset(doms[k], 3) for k in others)...))
    signs = Dict{Int,Float64}()
    function ratio(v)
        fa = _partial_n(F, v, i, doms[i])
        fb = _partial_n(F, v, j, doms[j])
        for (k, f) in ((i, fa), (j, fb))
            isfinite(f) || throw(
                ArgumentError(
                    "the result is not finite in :$(names[k]) somewhere in its domain",
                ),
            )
            iszero(f) && throw(
                ArgumentError(
                    "the result stops changing with :$(names[k]) somewhere in its " *
                    "domain, so it is not strictly monotone there",
                ),
            )
            s = get!(signs, k, sign(f))
            s == sign(f) || throw(
                ArgumentError("the result is not monotone in :$(names[k]) over its domain"),
            )
        end
        return log(abs(fa / fb))
    end
    worst = 0.0
    for a in inset(doms[i], n), b in inset(doms[j], n)
        v = copy(ref)
        v[i], v[j] = a, b
        ρ0 = ratio(v)
        for rest in restvalues
            w = copy(v)
            for (q, k) in enumerate(others)
                w[k] = rest[q]
            end
            worst = max(worst, abs(ratio(w) - ρ0))
        end
    end
    return worst
end

# Make inputs i and j the first stage, with the rest charted from their combination.
# Throws an ArgumentError saying, in a line, why this pair will not do.
function _pivot_on(
    F,
    vs::Vector{Variable},
    vz::Variable,
    i::Int,
    j::Int;
    width,
    height,
    depth,
)
    m = length(vs)
    doms = Tuple{Float64,Float64}[v.domain for v in vs]
    names = [v.name for v in vs]
    d = _pivot_defect(F, doms, names, i, j)
    d <= 1e-6 || throw(
        ArgumentError(
            "the result does not depend on them through a single combination of the " *
            "two (defect $(round(d; sigdigits = 3)))",
        ),
    )

    others = [k for k = 1:m if k != i && k != j]
    # Where each argument of F comes from: 0 for the pair, else its place among the rest.
    slot = Tuple(something(findfirst(==(k), others), 0) for k = 1:m)
    # The combination itself: F with every other input held at the middle of its domain.
    # Any fixed setting would do, since by the test above they all give the same level
    # curves; the scale it is drawn on is blank, so its units never show.
    ref = Tuple((doms[k][1] + doms[k][2]) / 2 for k in others)
    T = _combination(F, Val(m), i, j, slot, ref)
    tname = Symbol(:_t, depth)
    blank = Variable(tname; label = "", ticks = Float64[])
    # The first stage, drawn with either input on the left. Which goes where, and how wide
    # the stage is, are layout choices settled once the chain is assembled.
    function first_stage(w, swapped)
        swapped ||
            return _auto_two(T, vs[i], vs[j], blank; width = w, height, verify = false)
        return _auto_two(
            (b, a) -> T(a, b),
            vs[j],
            vs[i],
            blank;
            width = w,
            height,
            verify = false,
        )
    end
    stage = try
        first_stage(width, false)
    catch err
        err isa ArgumentError || rethrow()
        throw(
            ArgumentError(
                "their combination has no three-scale chart of its own, so it cannot be " *
                "the first stage",
            ),
        )
    end

    # The rest of the chain needs F as a function of the combination, G(t, rest…): some
    # point of the (xᵢ, xⱼ) rectangle where T equals t. Walk a path along which T rises
    # and tabulate the path parameter against t once, so each evaluation of G costs an
    # interpolation rather than a root find -- which matters when chains nest.
    P = _Partials(T, doms[i], doms[j])
    path = _rising_path(T, P, doms[i], doms[j])
    dsdt = [
        1 / (_Fx(P, path.along(s)...) * path.dx + _Fy(P, path.along(s)...) * path.dy)
        for s in path.svals
    ]
    s_of_t = TabulatedScale(path.zs, path.svals, dsdt)
    G = _remainder(F, Val(m), i, j, slot, path.along, s_of_t)

    tvar = Variable(tname; domain = stage[tname].domain, label = "", ticks = Float64[])
    downstream = try
        _chain_stages(G, Variable[tvar; vs[others]], vz; width, height, depth = depth + 1)
    catch err
        err isa ArgumentError || rethrow()
        throw(ArgumentError("what remains once they are combined cannot be charted"))
    end
    assembled = try
        compound(downstream..., stage)
    catch err
        err isa ArgumentError || rethrow()
        throw(
            ArgumentError(
                "the two stages want their combination on different scales, so this " *
                "pair needs a grid rather than a pivot",
            ),
        )
    end
    fitted = _fit_first_stage(first_stage, assembled, stage, downstream, tname, width)
    return Nomogram[downstream; fitted]
end

# F's arguments, with the pair at i and j and the rest in order elsewhere. A tuple of known
# length, so F is called with a fixed number of arguments rather than a splatted vector --
# which matters because the user's F sits behind `_Formula` and is called dynamically.
@inline _arguments(::Val{M}, i, j, slot, a, b, rest) where {M} =
    ntuple(k -> k == i ? a : k == j ? b : Float64(rest[slot[k]]), Val(M))

# The combination of the pair, and F as a function of it and the rest. Built behind these
# function barriers so that the arity is a type parameter.
_combination(F, v::Val, i, j, slot, ref) =
    (a, b) -> float(F(_arguments(v, i, j, slot, a, b, ref)...))

function _remainder(F, v::Val, i, j, slot, along, s_of_t)
    return function (t, rest...)
        a, b = along(clamp(transform(s_of_t, t), 0.0, 1.0))
        return float(F(_arguments(v, i, j, slot, a, b, rest)...))
    end
end

# The pivot line sits between the first stage's inputs, so one of them always lands
# inside the rest of the chain -- and left to chance, on top of one of its axes. So choose
# which input goes inside and how wide the stage is, aiming the inside input at the middle
# of a gap between the axes already there, and keep whichever choice leaves the new axes
# furthest from their neighbours for the width it takes -- otherwise a lopsided stage
# would fling its outside input a long way to buy a little room inside. An input that is
# itself a pivot line goes outside: the stage feeding it straddles it in turn and needs
# the room. The stage is never wider than the
# last one, nor narrower than a quarter of it. The pivot line always reads upwards, so
# whichever way the stages are aligned, swapping the inputs mirrors the stage about it.
function _fit_first_stage(first_stage, assembled, stage, downstream, tname, width)
    consumer = only(st for st in downstream if tname in variablenames(st))
    x(n) = (s = assembled[n]; point(s, sum(s.domain) / 2)[1])
    xt = x(tname)
    inward = sign(x(name(consumer.scales[3])) - xt)
    # Everything as signed offsets from the pivot line, positive into the chain.
    # A pivot line belongs to two stages; it is one axis.
    others = unique(name(s) for st in downstream for s in st.scales)
    placed = sort!([inward * (x(n) - xt) for n in others])
    inputs = [name(s) for s in stage.scales[1:2]]
    offsets = [inward * (x(n) - xt) for n in inputs]
    pivots = [_is_pivot_line(n) for n in inputs]
    gaps = [(placed[k] + placed[k+1]) / 2 for k = 1:(length(placed)-1) if placed[k] >= 0]

    best = (spread = -Inf, w = stage.width, swapped = false)
    for swapped in (false, true)
        o = swapped ? -offsets : offsets
        inside = findfirst(>(0), o)
        isnothing(inside) && continue
        # Mirror so that a pivot-line input stays outside, unless both are.
        any(pivots) && !all(pivots) && pivots[inside] && continue
        for target in gaps
            w = clamp(stage.width * target / o[inside], width / 4, width)
            # Only the gaps beside the new axes: the rest are already settled.
            new = o .* (w / stage.width)
            spread =
                min(abs(new[1] - new[2]), minimum(abs(p - q) for p in new, q in placed))
            spread /=
                max(maximum(placed), maximum(new)) - min(minimum(placed), minimum(new))
            if spread > best.spread + 1e-9
                best = (; spread, w, swapped)
            end
        end
    end
    (best.swapped || !isapprox(best.w, stage.width; rtol = 1e-3)) || return stage
    return first_stage(best.w, best.swapped)
end

_is_pivot_line(n::Symbol) = occursin(r"^_t\d+$", string(n))

# The narrowest gap between neighbouring axes of a chain, as a fraction of its width.
function _spacing(stages)
    nom = compound(stages...)
    xs = sort!([point(s, sum(s.domain) / 2)[1] for s in nom.scales])
    return minimum(diff(xs)) / (xs[end] - xs[1])
end

# The stages of a chain for z = F(vs…), in an order compound accepts: every stage shares a
# variable with one placed before it. Each pivot is built after the stages that consume
# it, so reversing the order of construction gives exactly that.
#
# Every pairing that works is built and the best spaced kept. Which pair comes first
# decides which ranges share the last stage, and a result axis squeezed between a scale
# spanning decades and one spanning a fraction of one sits almost on the former; no
# choice of widths can fix that afterwards. Below the second level the first pairing that
# works is taken: the search multiplies with every level, and the first two levels settle
# both the last stage and whether the chain branches.
function _chain_stages(F, vs::Vector{Variable}, vz::Variable; width, height, depth)
    length(vs) == 2 &&
        return Nomogram[_auto_two(F, vs[1], vs[2], vz; width, height, verify = false)]
    reasons = String[]
    best, bestspacing = nothing, -Inf
    for i = 1:(length(vs)-1), j = (i+1):length(vs)
        try
            stages = _pivot_on(F, vs, vz, i, j; width, height, depth)
            spacing = _spacing(stages)
            if spacing > bestspacing + 1e-9
                best, bestspacing = stages, spacing
            end
            depth > 2 && break
        catch err
            err isa ArgumentError || rethrow()
            push!(reasons, ":$(vs[i].name) with :$(vs[j].name) — $(err.msg)")
        end
    end
    isnothing(best) || return best
    throw(
        ArgumentError(
            "cannot split $(vz.name) = F($(join(string.(getfield.(vs, :name)), ", "))) " *
            "into a chain of three-scale charts. Each pair of inputs was tried as the " *
            "first stage:\n" *
            join(("  • " * r for r in reasons), "\n"),
        ),
    )
end

function _auto_chain(F, vs::Vector{Variable}, vz::Variable; width, height, title, rtol)
    stages = _chain_stages(F, vs, vz; width, height, depth = 1)
    nom = compound(stages...; title)

    # As with two inputs, the chart is only worth returning if reading it gives F.
    zd = nom[vz.name].domain
    zspan = zd[2] - zd[1]
    names = Tuple(v.name for v in vs)
    worst = 0.0
    for vals in Iterators.product((range(v.domain...; length = 4) for v in vs)...)
        got = _evaluate(nom; NamedTuple{names}(vals)...)[vz.name]
        worst = max(worst, abs(got - float(F(vals...))) / zspan)
    end
    worst <= rtol || throw(
        ArgumentError(
            "the derived chart disagrees with F by $(round(worst; sigdigits = 3)) of the " *
            "result's range, more than rtol = $rtol, so it has not been returned. This " *
            "should not happen for a formula that splits into separable stages; please " *
            "report it.",
        ),
    )
    return nom
end
