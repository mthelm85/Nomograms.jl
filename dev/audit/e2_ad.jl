# Experiment 2: finite differences vs ForwardDiff
using Pkg;
Pkg.activate(joinpath(@__DIR__, "env"); io = devnull)
using Nomograms, ForwardDiff, Printf, Unitful
fmt(x) = @sprintf("%.2e", x)

function defect_ad(F, xd, yd; n = 9)
    inset(d) = range(d[1] + 0.02 * (d[2] - d[1]), d[2] - 0.02 * (d[2] - d[1]); length = n)
    xs, ys = inset(float.(xd)), inset(float.(yd))
    L = [
        log(
            abs(
                ForwardDiff.derivative(t -> F(t, y), x) /
                ForwardDiff.derivative(t -> F(x, t), y),
            ),
        ) for x in xs, y in ys
    ]
    r = (n + 1) ÷ 2
    maximum(abs(L[i, j] - L[i, r] - L[r, j] + L[r, r]) for i = 1:n, j = 1:n)
end

println("== noise floor on separable formulas: FD vs AD ==")
for (n, F, xd, yd) in [
    ("x*y", (x, y) -> x * y, (1, 10), (1, 10)),
    ("x+y", (x, y) -> x + y, (0, 10), (0, 10)),
    ("hypot", (x, y) -> hypot(x, y), (1, 10), (1, 10)),
    ("x*y wide", (x, y) -> x * y, (1e-3, 1e3), (1e-3, 1e3)),
    ("exp(10x)*y", (x, y) -> exp(10x) * y, (0, 2), (1, 10)),
    ("1e-9*(x+y)", (x, y) -> 1e-9 * (x + y), (0, 10), (0, 10)),
]
    println(
        rpad(n, 12),
        " FD ",
        fmt(separability_defect(F, xd, yd)),
        "   AD ",
        fmt(defect_ad(F, xd, yd)),
    )
end

println("\n== recovered transform vs exact, f(x) = log x for x*y ==")
# Rebuild the package's integration with each derivative and compare to log up to an affine map.
function recovered(deriv, lo, hi)
    t = Nomograms._integrate_to_scale(deriv, lo, hi)
    xs = exp.(range(log(lo), log(hi); length = 400))
    got = [Nomograms.transform(t, x) for x in xs]
    want = log.(xs) .- log(lo)
    # fit got ≈ a*want (both start at 0)
    a = sum(got .* want) / sum(want .^ 2)
    maximum(abs.(got .- a .* want)) / (want[end] - want[1])
end
for (lo, hi) in ((1.0, 10.0), (1.0, 1000.0), (1e-3, 1e3))
    F(x, y) = x * y
    yr = sqrt(lo * hi)
    # The package's f'(x): F_x / F_y along y = yr, each at (x, yr).
    fd(x) =
        Nomograms._deriv(t -> F(t, yr), x, lo, hi) /
        Nomograms._deriv(t -> F(x, t), yr, lo, hi)
    ad(x) =
        ForwardDiff.derivative(t -> F(t, yr), x) / ForwardDiff.derivative(t -> F(x, t), yr)
    println(
        "domain [",
        lo,
        ", ",
        hi,
        "]   FD rel.err ",
        fmt(recovered(fd, lo, hi)),
        "   AD rel.err ",
        fmt(recovered(ad, lo, hi)),
    )
end

println("\n== can AD run through the package's formula wrapper and units? ==")
w = Nomograms._Formula((x, y) -> x * y, 2)
try
    ForwardDiff.derivative(t -> w(t, 2.0), 3.0)
    println("wrapper: ok")
catch e
    println("wrapper: fails -- ", first(sprint(showerror, e), 110))
end
try
    g(t) = ustrip(u"kg/m^2", (t * u"kg") / (1.75u"m")^2)
    println("unitful: ", ForwardDiff.derivative(g, 70.0))
catch e
    println("unitful: fails -- ", first(sprint(showerror, e), 110))
end

println("\n== cost: time per defect evaluation ==")
F(x, y) = x * y
separability_defect(F, (1, 10), (1, 10));
defect_ad(F, (1, 10), (1, 10))
t_fd = @elapsed for _ = 1:200
    separability_defect(F, (1, 10), (1, 10))
end
t_ad = @elapsed for _ = 1:200
    defect_ad(F, (1, 10), (1, 10))
end
println("FD ", fmt(t_fd / 200), " s   AD ", fmt(t_ad / 200), " s")
