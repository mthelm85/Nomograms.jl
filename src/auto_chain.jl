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
    function arguments(a, b, rest)
        v = Vector{Float64}(undef, m)
        v[i], v[j] = a, b
        for (q, k) in enumerate(others)
            v[k] = rest[q]
        end
        return v
    end
    # The combination itself: F with every other input held at the middle of its domain.
    # Any fixed setting would do, since by the test above they all give the same level
    # curves; the scale it is drawn on is blank, so its units never show.
    ref = [(doms[k][1] + doms[k][2]) / 2 for k in others]
    T = (a, b) -> float(F(arguments(a, b, ref)...))
    tname = Symbol(:_t, depth)
    blank = Variable(tname; label = "", ticks = Float64[])
    # The first stage, drawn with either input on the left. Which goes where, and how wide
    # the stage is, are layout choices settled once the chain is assembled.
    function first_stage(w, swapped)
        swapped || return _auto_two(T, vs[i], vs[j], blank; width = w, height)
        return _auto_two((b, a) -> T(a, b), vs[j], vs[i], blank; width = w, height)
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
    G = function (t, rest...)
        a, b = path.along(clamp(transform(s_of_t, t), 0.0, 1.0))
        return float(F(arguments(a, b, rest)...))
    end

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
                "pair needs a net rather than a pivot",
            ),
        )
    end
    fitted = _fit_first_stage(first_stage, assembled, stage, downstream, tname, width)
    return Nomogram[downstream; fitted]
end

# The turning line sits between the first stage's inputs, so one of them always lands
# inside the stage it feeds -- and left to chance, on top of that stage's result axis. So
# put the input nearer the turning line on that side, and size the stage so it lands
# halfway to the result, never wider than the stage it feeds. The turning line always
# reads upwards, so whichever way the stages are aligned, swapping the inputs mirrors the
# stage about it.
function _fit_first_stage(first_stage, assembled, stage, downstream, tname, width)
    consumer = only(st for st in downstream if tname in variablenames(st))
    x(n) = (s = assembled[n]; point(s, sum(s.domain) / 2)[1])
    xt = x(tname)
    xr = x(name(consumer.scales[3]))
    # Signed offsets from the turning line, positive towards the next stage's result.
    da, db = (sign(xr - xt) * (x(name(s)) - xt) for s in stage.scales[1:2])
    near = min(abs(da), abs(db))
    swapped = max(da, db) > near
    w = min(stage.width * abs(xr - xt) / 2 / near, width)
    (swapped || !isapprox(w, stage.width; rtol = 1e-3)) || return stage
    return first_stage(w, swapped)
end

# The stages of a chain for z = F(vs…), in an order compound accepts: every stage shares a
# variable with one placed before it. Each pivot is built after the stages that consume
# it, so reversing the order of construction gives exactly that.
function _chain_stages(F, vs::Vector{Variable}, vz::Variable; width, height, depth)
    length(vs) == 2 && return Nomogram[_auto_two(F, vs[1], vs[2], vz; width, height)]
    reasons = String[]
    for i = 1:(length(vs)-1), j = (i+1):length(vs)
        try
            return _pivot_on(F, vs, vz, i, j; width, height, depth)
        catch err
            err isa ArgumentError || rethrow()
            push!(reasons, ":$(vs[i].name) with :$(vs[j].name) — $(err.msg)")
        end
    end
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
        got = evaluate(nom; NamedTuple{names}(vals)...)[vz.name]
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
