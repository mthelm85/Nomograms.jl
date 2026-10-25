using Pkg;
Pkg.activate(joinpath(@__DIR__, "env"); io = devnull)
using Nomograms, Printf
fmt(x) = @sprintf("%.2e", x)
pd(F, doms, i, j) = Nomograms._pivot_defect(
    F,
    Tuple{Float64,Float64}[Float64.(d) for d in doms],
    [:x, :y, :w][1:length(doms)],
    i,
    j,
)
println("== pairs that do combine (should be ~0) ==")
for (n, F, doms, i, j) in [
    ("x*y*w (x,y)", (x, y, w) -> x * y * w, [(1, 10), (1, 10), (1, 10)], 1, 2),
    ("(x+y)*w (x,y)", (x, y, w) -> (x + y) * w, [(1, 10), (1, 10), (1, 10)], 1, 2),
    (
        "n R T / V (n,T)",
        (n, T, V) -> n * 8.314 * T / V,
        [(0.5, 5), (250, 400), (0.01, 0.1)],
        1,
        2,
    ),
    ("exp(x+y)*w", (x, y, w) -> exp(x + y) * w, [(0, 3), (0, 3), (1, 10)], 1, 2),
    ("hypot(x,y)+w", (x, y, w) -> hypot(x, y) + w, [(1, 10), (1, 10), (1, 10)], 1, 2),
]
    println(rpad(n, 18), fmt(pd(F, doms, i, j)))
end
println("\n== pairs that do not (should be large) ==")
for (n, F, doms, i, j) in [
    ("(x+y)*w (x,w)", (x, y, w) -> (x + y) * w, [(1, 10), (1, 10), (1, 10)], 1, 3),
    ("x*y+w*x (x,y)", (x, y, w) -> x * y + w * x^2, [(1, 10), (1, 10), (1, 10)], 1, 2),
]
    println(rpad(n, 18), fmt(pd(F, doms, i, j)))
end
println("\n== x*y*w + eps*x*w^2, pair (x,y) ==")
for e in (1e-2, 1e-4, 1e-6, 1e-7, 1e-8)
    println(
        "eps=",
        fmt(e),
        "  defect=",
        fmt(pd((x, y, w) -> x * y * w + e * x * w^2, [(1, 10), (1, 10), (1, 10)], 1, 2)),
    )
end
