"""
    CompoundNomogram

Several [`Nomogram`](@ref) stages sharing variables, for relations that need more than
three variables.

Three scales is the limit of a single alignment line, so a four-variable relation is split
into two: `g(a, b) = t` and `h(t, c) = d`, sharing the intermediate `t`. The reader lays a
line across the first stage to find where it crosses the shared scale, then pivots about
that point and reads the second — which is why the shared scale is often called the pivot
or turning line.

Build one with [`compound`](@ref), which positions the stages so their shared scales
coincide. Drawing works exactly as for a single nomogram; solving propagates values from
stage to stage (see [`evaluate`](@ref)).
"""
struct CompoundNomogram <: AbstractNomogram
    stages::Vector{Nomogram}
    scales::Vector{Scale}
    width::Float64
    height::Float64
    title::String
end

# The similarity (rotate, scale, translate) carrying one segment onto another. Complex
# arithmetic does the rotate-and-scale in a single division. A similarity is affine, so it
# cannot disturb collinearity and therefore cannot change what a stage computes.
function _similarity(src_a::Point, src_b::Point, dst_a::Point, dst_b::Point)
    z_src = complex(src_b[1] - src_a[1], src_b[2] - src_a[2])
    z_dst = complex(dst_b[1] - dst_a[1], dst_b[2] - dst_a[2])
    abs(z_src) > 1e-12 ||
        throw(ArgumentError("cannot align on a shared scale of zero length"))
    s = z_dst / z_src
    return function (p)
        z = s * complex(p[1] - src_a[1], p[2] - src_a[2])
        return (real(z) + dst_a[1], imag(z) + dst_a[2])
    end
end

function _translated(s::Scale, f)
    curve = s.curve
    return Scale(s.variable, s.domain, t -> f(Point(curve(float(t)))), s.side)
end

# Line up `stage` on an already-placed scale for the variable they share.
function _align_stage(stage::Nomogram, placed::Scale, shared::Symbol, atol::Float64)
    mine = scale(stage, shared)
    lo, hi = mine.domain
    plo, phi = placed.domain
    if !(
        isapprox(lo, plo; rtol = 1e-6, atol = 1e-12) &&
        isapprox(hi, phi; rtol = 1e-6, atol = 1e-12)
    )
        throw(
            ArgumentError(
                "the stages disagree about :$shared: one spans " *
                "($(default_tickformat(plo)), $(default_tickformat(phi))) and the " *
                "other ($(default_tickformat(lo)), $(default_tickformat(hi))). A " *
                "shared scale has to carry the same values in both, so give the " *
                "second stage `domain = ($(default_tickformat(plo)), " *
                "$(default_tickformat(phi)))` for :$shared -- the range the first " *
                "stage derived.",
            ),
        )
    end
    a, b = lo + (hi - lo) / 4, lo + 3 * (hi - lo) / 4
    f = _similarity(point(mine, a), point(mine, b), point(placed, a), point(placed, b))

    # Matching two values lines up a straight scale completely, but nothing guarantees the
    # shared scale is straight or parametrised alike in both stages. Check the rest.
    span = max(hypot((point(placed, hi) .- point(placed, lo))...), 1.0)
    for t in range(lo, hi; length = 17)
        got = f(point(mine, t))
        want = point(placed, t)
        if hypot(got[1] - want[1], got[2] - want[2]) > atol * span
            throw(
                ArgumentError(
                    "the two stages' scales for :$shared cannot be laid on top of " *
                    "one another: at $(default_tickformat(t)) they stay " *
                    "$(round(hypot(got[1] - want[1], got[2] - want[2]); digits = 3)) " *
                    "apart after alignment. They must be the same shape, which for " *
                    "the usual constructions means the same transform as well as the " *
                    "same domain.",
                ),
            )
        end
    end
    return f
end

"""
    compound(stages::Nomogram...; title = "", atol = 1e-6) -> CompoundNomogram

Join nomogram stages into one chart, positioning each so that the scales it shares with an
earlier stage lie exactly on top of them.

Variables are shared by name. Each stage after the first must share at least one with a
stage already placed, and that shared scale must carry the same values over the same
domain in both — so build the later stage with the domain the earlier one derived:

```julia
first = additive_nomogram(a, b, Variable(:t))          # :t derived
second = additive_nomogram(Variable(:t; domain = first[:t].domain), c, d)
chart = compound(first, second)
```

Alignment is a similarity — rotate, scale, translate — chosen to carry the later stage's
shared scale onto the placed one, then checked along the whole domain rather than just at
the two points that fixed it. A similarity is affine, so it cannot change what any stage
computes.

The shared scale is drawn once. Note that the stages generally interleave rather than
sitting side by side: a parallel-axis construction puts its result *between* its two
inputs, so a pivot shared by two such stages has both stages' other scales straddling it.
[`project`](@ref) will reshape the result if the arrangement is awkward.
"""
function compound(stages::Nomogram...; title::AbstractString = "", atol::Real = 1e-6)
    length(stages) >= 2 ||
        throw(ArgumentError("a compound nomogram needs at least 2 stages"))
    placed = Nomogram[stages[1]]
    byname = Dict{Symbol,Scale}(name(s) => s for s in stages[1].scales)
    order = Scale[stages[1].scales...]

    for stage in stages[2:end]
        shared = [n for n in variablenames(stage) if haskey(byname, n)]
        isempty(shared) && throw(
            ArgumentError(
                "stage $(length(placed) + 1) shares no variable with the stages " *
                "before it (it has " *
                join(":" .* string.(variablenames(stage)), ", ") *
                "); the stages of a compound nomogram are linked by a shared scale",
            ),
        )
        f = _align_stage(stage, byname[shared[1]], shared[1], float(atol))
        moved = Scale[_translated(s, f) for s in stage.scales]
        push!(
            placed,
            Nomogram(
                moved;
                width = stage.width,
                height = stage.height,
                title = stage.title,
                relation = stage.relation,
            ),
        )
        for s in moved
            haskey(byname, name(s)) && continue
            byname[name(s)] = s
            push!(order, s)
        end
    end

    # Put the union of everything drawn at the origin, so the compound chart has a canvas
    # of its own rather than inheriting one stage's.
    xs = Float64[]
    ys = Float64[]
    for s in order, t in range(s.domain...; length = 65)
        p = point(s, t)
        push!(xs, p[1])
        push!(ys, p[2])
    end
    x0, y0 = minimum(xs), minimum(ys)
    shift(p) = (p[1] - x0, p[2] - y0)
    # Left to right, so the list of scales reads in the order they appear on the page.
    # Stage order would interleave, because a parallel construction puts its result
    # between its inputs and the stages therefore straddle the shared scale.
    sort!(order; by = s -> (point(s, (s.domain[1] + s.domain[2]) / 2)[1], string(name(s))))
    final = Scale[_translated(s, shift) for s in order]
    stagesout = Nomogram[
        Nomogram(
            Scale[_translated(s, shift) for s in st.scales];
            width = st.width,
            height = st.height,
            title = st.title,
            relation = st.relation,
        ) for st in placed
    ]
    return CompoundNomogram(
        stagesout,
        final,
        maximum(xs) - x0,
        maximum(ys) - y0,
        String(title),
    )
end

"""
    evaluate(cn::CompoundNomogram; values...) -> NamedTuple

Solve a compound nomogram by reading its stages in turn.

Supply enough variables that some stage has exactly one unknown; that stage is solved, its
answer becomes known, and the pass repeats until everything is determined. Which variables
you supply is up to you, as long as the chain closes — the intermediate on a shared scale
is usually left to be derived, but it can equally be an input.

Throws, naming what is left, if the values given are not enough to finish.
"""
function evaluate(cn::CompoundNomogram; values...)
    names = variablenames(cn)
    given = Dict{Symbol,Float64}()
    for (k, v) in pairs(values)
        k in names || throw(
            ArgumentError(
                "no variable :$k in this nomogram (have " *
                join(":" .* string.(names), ", ") *
                ")",
            ),
        )
        given[k] = float(v)
    end
    all(n -> haskey(given, n), names) &&
        throw(ArgumentError("every variable was given, so there is nothing to solve"))

    progress = true
    while progress
        progress = false
        for stage in cn.stages
            sn = variablenames(stage)
            missingnames = [n for n in sn if !haskey(given, n)]
            length(missingnames) == 1 || continue
            knownnames = Tuple(n for n in sn if n != missingnames[1])
            args = NamedTuple{knownnames}(Tuple(given[n] for n in knownnames))
            sol = evaluate(stage; args...)
            given[missingnames[1]] = sol[missingnames[1]]
            progress = true
        end
    end

    left = [n for n in names if !haskey(given, n)]
    isempty(left) || throw(
        ArgumentError(
            "cannot finish from the values given: " *
            join(":" .* string.(left), ", ") *
            " stayed unknown. Every stage needs all but one of its variables before " *
            "it can be read, so supply another input.",
        ),
    )
    return NamedTuple{Tuple(names)}(Tuple(given[n] for n in names))
end

"""
    isopleth(cn::CompoundNomogram; label = "", values...) -> Vector{Isopleth}
    isopleth(cn::CompoundNomogram, solution::NamedTuple; label = "") -> Vector{Isopleth}

The alignment lines for one reading of a compound nomogram: one per stage, since reading
it means drawing a line across each in turn.

They come back sharing a `group`, so a renderer draws the whole reading in one colour. Pass
the vector straight to [`drawing`](@ref) or a backend. The label goes on the last stage,
where the answer is read.
"""
function isopleth(
    cn::CompoundNomogram,
    solution::NamedTuple;
    label::AbstractString = "",
    group::Int = 1,
)
    out = Isopleth[]
    for (i, stage) in enumerate(cn.stages)
        sol = NamedTuple{Tuple(variablenames(stage))}(
            Tuple(solution[n] for n in variablenames(stage)),
        )
        text = i == length(cn.stages) ? String(label) : ""
        push!(out, isopleth(stage, sol; label = text, group))
    end
    return out
end

function isopleth(
    cn::CompoundNomogram;
    label::AbstractString = "",
    group::Int = 1,
    values...,
)
    return isopleth(cn, evaluate(cn; values...); label, group)
end

function Base.show(io::IO, cn::CompoundNomogram)
    print(io, "CompoundNomogram(")
    join(io, (":$(name(s))" for s in cn.scales), ", ")
    print(io, "; ", length(cn.stages), " stages)")
    return nothing
end
