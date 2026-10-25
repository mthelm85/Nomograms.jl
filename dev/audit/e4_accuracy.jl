# Experiment 4: is the built-in check (9x9 / 4^n grid, error relative to the result's range)
# a good guard? Compare it with dense sampling under three measures.
using Pkg;
Pkg.activate(joinpath(@__DIR__, "env"); io = devnull)
using Nomograms, Printf, Random
fmt(x) = @sprintf("%.1e", x)
Random.seed!(1)

function measure(F, doms, names, zname)
    args = [n => d for (n, d) in zip(names, doms)]
    nom = auto_nomogram(F, args..., zname)
    zs = nom[zname]
    zd = zs.domain
    span = zd[2] - zd[1]
    # page position of a value on the result axis
    pos(z) = Nomograms.point(zs, z)
    axislen = hypot((pos(zd[2]) .- pos(zd[1]))...)
    # the package's own grid
    k = length(doms) == 2 ? 9 : 4
    grid = Iterators.product((range(d...; length = k) for d in doms)...)
    pkg = maximum(
        abs(Nomograms._evaluate(nom; NamedTuple{Tuple(names)}(v)...)[zname] - F(v...)) /
        span for v in grid
    )
    # dense random sampling, including near the domain edges
    worst_range = worst_rel = worst_page = 0.0
    for _ = 1:3000
        v = Tuple(d[1] + rand()^(rand() < 0.3 ? 3 : 1) * (d[2] - d[1]) for d in doms)
        got = Nomograms._evaluate(nom; NamedTuple{Tuple(names)}(v)...)[zname]
        want = F(v...)
        worst_range = max(worst_range, abs(got - want) / span)
        worst_rel = max(worst_rel, abs(got / want - 1))
        worst_page =
            max(worst_page, hypot((pos(got) .- pos(clamp(want, zd...))) ...) / axislen)
    end
    return pkg, worst_range, worst_rel, worst_page
end

cases = [
    ("x*y", (x, y) -> x * y, [(1, 10), (1, 10)]),
    ("BMI", (m, h) -> m / h^2, [(40, 150), (1.4, 2.1)]),
    ("x*y 6 decades", (x, y) -> x * y, [(1e-3, 1e3), (1e-3, 1e3)]),
    ("exp(10x)*y", (x, y) -> exp(10x) * y, [(0, 2), (1, 10)]),
    ("x^10*y", (x, y) -> x^10 * y, [(1, 3), (1, 10)]),
    ("(1+r)^n", (r, n) -> (1 + r / 100)^n, [(1, 20), (1, 30)]),
    ("atan(y,x)", (y, x) -> atan(y, x), [(1, 10), (1, 10)]),
    ("parallel R", (a, b) -> 1 / (1 / a + 1 / b), [(1, 10), (1, 10)]),
    ("gas (3 inputs)", (n, T, V) -> n * 8.314 * T / V, [(0.5, 5), (250, 400), (0.01, 0.1)]),
    ("product of 4", (a, b, c, d) -> a * b * c * d, [(1, 5), (1, 5), (1, 5), (1, 5)]),
]
println(
    rpad("formula", 16),
    rpad("pkg check", 11),
    rpad("dense/range", 13),
    rpad("dense rel", 12),
    "dense page",
)
for (n, F, doms) in cases
    names = [Symbol(:v, i) for i = 1:length(doms)]
    try
        pkg, wr, wrel, wp = measure(F, doms, names, :z)
        println(
            rpad(n, 16),
            rpad(fmt(pkg), 11),
            rpad(fmt(wr), 13),
            rpad(fmt(wrel), 12),
            fmt(wp),
        )
    catch e
        println(rpad(n, 16), "refused: ", first(sprint(showerror, e), 80))
    end
end
