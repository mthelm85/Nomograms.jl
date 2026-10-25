# Experiment 5: does root finding miss crossings? Adversarial readings on the curved scale.
using Pkg;
Pkg.activate(joinpath(@__DIR__, "env"); io = devnull)
using Nomograms, Printf
fmt(x) = @sprintf("%.1e", x)
quad = determinant_nomogram(
    Variable(:p; domain = (-12, 0)) => p -> (0.0, p),
    Variable(:q; domain = (0, 36)) => q -> (1.0, q),
    Variable(:z; domain = (0.2, 6)) => z -> (1 / (1 + z), -z^2 / (1 + z));
    relation = (p, q, z) -> z^2 + p * z + q,
)
found(p, q) =
    try
        sort([s.z for s in evaluate_all(quad; p, q)])
    catch
        Float64[]
    end

println("== two roots r and r + d, r = 2 ==")
for d in (1.0, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8, 1e-10)
    r1, r2 = 2.0, 2.0 + d
    z = found(-(r1 + r2), r1 * r2)
    err = length(z) == 2 ? max(abs(z[1] - r1), abs(z[2] - r2)) : NaN
    println(
        "d=",
        rpad(fmt(d), 9),
        " found ",
        length(z),
        length(z) == 2 ? "  max err " * fmt(err) : "   <-- expected 2",
    )
end
println("\n== a double root, and lines just missing the curve (should find 0) ==")
for e in (0.0, 1e-12, 1e-8, 1e-4, 1e-2)
    z = found(-4.0, 4.0 + e)      # z^2 - 4z + (4+e): double root at 2 when e = 0
    println(
        "q = 4 + ",
        rpad(fmt(e), 8),
        " found ",
        length(z),
        length(z) > 0 ? "  at " * join(fmt.(z), ", ") : "",
    )
end
println("\n== roots near the ends of the scale ==")
for (r1, r2) in ((0.2, 3.0), (0.2000001, 3.0), (0.21, 0.22), (5.9, 6.0), (5.99999, 6.0))
    z = found(-(r1 + r2), r1 * r2)
    println("roots ", r1, ", ", r2, "  found ", length(z), "  ", join(fmt.(z), ", "))
end
println("\n== scan: random pairs of roots inside the scale, how many are missed? ==")
using Random;
Random.seed!(2)
miss = 0;
total = 0;
closest_missed = Inf
for _ = 1:4000
    r1 = 0.2 + 5.8 * rand()
    r2 = 0.2 + 5.8 * rand()
    abs(r1 - r2) < 1e-12 && continue
    p, q = -(r1 + r2), r1 * r2
    (-12 <= p <= 0 && 0 <= q <= 36) || continue
    global total += 1
    if length(found(p, q)) != 2
        global miss += 1
        global closest_missed = min(closest_missed, abs(r1 - r2))
    end
end
println(
    "missed ",
    miss,
    " of ",
    total,
    "  (closest gap among misses: ",
    fmt(closest_missed),
    ")",
)
println("\ntime per evaluate_all: ", fmt(@elapsed(for _ = 1:500
    found(-5.0, 6.0)
end) / 500), " s")
