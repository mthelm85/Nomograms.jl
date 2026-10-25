# Experiment 1: where does separability_defect fall, and does 1e-6 separate the classes?
using Pkg;
Pkg.activate(joinpath(@__DIR__, "env"); io = devnull)
using Nomograms, Printf
sd(F, xd, yd) = separability_defect(F, xd, yd)
fmt(x) = @sprintf("%.2e", x)

println("== separable formulas (should be ~0) ==")
sep = [
    ("x*y", (x, y) -> x * y, (1, 10), (1, 10)),
    ("x/y", (x, y) -> x / y, (1, 10), (1, 10)),
    ("x+y", (x, y) -> x + y, (0, 10), (0, 10)),
    ("x^2*y^3", (x, y) -> x^2 * y^3, (1, 10), (1, 10)),
    ("exp(x)*y", (x, y) -> exp(x) * y, (0, 5), (1, 10)),
    ("(1+x)(1+y)-1", (x, y) -> x + y + x * y, (0, 5), (0, 5)),
    ("hypot", (x, y) -> hypot(x, y), (1, 10), (1, 10)),
    ("atan(y,x)", (x, y) -> atan(y, x), (1, 10), (1, 10)),
    ("x^y", (x, y) -> x^y, (1.5, 10), (0.5, 3)),
    ("parallel R", (x, y) -> 1 / (1 / x + 1 / y), (1, 10), (1, 10)),
    ("BMI", (m, h) -> m / h^2, (40, 150), (1.4, 2.1)),
    ("x*y wide", (x, y) -> x * y, (1e-3, 1e3), (1e-3, 1e3)),
    ("x*y tiny dom", (x, y) -> x * y, (1e-7, 1e-6), (1e-7, 1e-6)),
    ("x*y huge dom", (x, y) -> x * y, (1e6, 1e7), (1e6, 1e7)),
    ("x+y near 0", (x, y) -> x + y, (0, 1e-3), (0, 1e-3)),
    ("1e6*x*y", (x, y) -> 1e6 * x * y, (1, 10), (1, 10)),
    ("x*y + 1e6", (x, y) -> x * y + 1e6, (1, 10), (1, 10)),
    ("1e-9*(x+y)", (x, y) -> 1e-9 * (x + y), (0, 10), (0, 10)),
    ("exp(10x)*y", (x, y) -> exp(10x) * y, (0, 2), (1, 10)),
    ("x^10*y", (x, y) -> x^10 * y, (1, 3), (1, 10)),
]
worst_sep = 0.0
for (n, F, xd, yd) in sep
    d = sd(F, xd, yd)
    global worst_sep = max(worst_sep, d)
    println(rpad(n, 16), fmt(d), d > 1e-6 ? "   <-- WRONGLY REFUSED" : "")
end
println("worst separable: ", fmt(worst_sep))

println("\n== x + y + eps*x*y^2 : not separable for eps != 0 ==")
for e in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8)
    d = sd((x, y) -> x + y + e * x * y^2, (0.5, 5), (0.5, 5))
    println(
        "eps=",
        fmt(e),
        "  defect=",
        fmt(d),
        d <= 1e-6 ? "   <-- accepted as separable" : "",
    )
end
println(
    "\n== x*y*(1 + eps*x) * ... : x*y + eps*x^2*y? (eps*x^2*y = x^2 y is separable!) use x*y + eps*(x+y)^2 ==",
)
for e in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8)
    d = sd((x, y) -> x * y + e * (x + y)^2, (1, 10), (1, 10))
    println(
        "eps=",
        fmt(e),
        "  defect=",
        fmt(d),
        d <= 1e-6 ? "   <-- accepted as separable" : "",
    )
end
println("\n== genuinely non-separable ==")
for (n, F, xd, yd) in [
    ("x^2+xy+y", (x, y) -> x^2 + x * y + y, (0.5, 5), (0.5, 5)),
    ("sin(x)+exp(y)+x*y", (x, y) -> sin(x) + exp(y) + x * y, (0.1, 1), (0.1, 1)),
    ("(x+y)^2+x", (x, y) -> (x + y)^2 + x, (1, 2), (1, 2)),
]
    println(rpad(n, 18), fmt(sd(F, xd, yd)))
end

# When a near-separable formula IS accepted, does the final verification catch the error?
println("\n== does verification catch an accepted near-separable formula? ==")
for e in (1e-6, 1e-7, 1e-8)
    try
        nom =
            auto_nomogram((x, y) -> x + y + e * x * y^2, :x => (0.5, 5), :y => (0.5, 5), :z)
        err = maximum(
            abs(evaluate(nom; x = a, y = b).z - (a + b + e * a * b^2)) for
            a in range(0.5, 5; length = 60), b in range(0.5, 5; length = 60)
        )
        span = -(reverse(nom[:z].domain)...)
        println(
            "eps=",
            fmt(e),
            "  chart built; worst reading error / range = ",
            fmt(err / span),
        )
    catch err
        println("eps=", fmt(e), "  refused: ", first(sprint(showerror, err), 90))
    end
end
