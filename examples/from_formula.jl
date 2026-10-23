# Charts from formulas as you would ordinarily write them.
#
# Every other example rewrites its formula first: takes logs, spots a product, picks the
# transforms, supplies their inverses, works out which variable goes in the middle.
# `auto_nomogram` does that part. Hand it the function and the ranges, and it finds the
# transforms -- or tells you there are none.

using Nomograms

# Compound growth, straight from the formula. Compare `compound_growth.jl`, which derives
# the same chart by hand: logs, a product, a custom transform with its inverse, LogScale(ℯ).
growth = auto_nomogram(
    (r, n) -> (1 + r / 100)^n,
    Variable(:rate; domain = (1, 20), label = "rate (% p.a.)"),
    Variable(:years; domain = (1, 30), label = "years"),
    # The growth axis compresses hard -- position goes as log(log(growth)) -- so automatic
    # ticks would label little above 10. Explicit ticks are kept, and the ones that do not
    # fit are still dropped.
    Variable(
        :growth;
        label = "growth multiple",
        ticks = [1.02, 1.05, 1.1, 1.25, 1.5, 2, 3, 5, 10, 20, 50, 100],
    );
    title = "Compound growth",
)
@show evaluate(growth; rate = 5, years = 10).growth      # 1.05^10
@show evaluate(growth; rate = 7, growth = 2).years       # the rule of 72's question

# Two resistors in parallel, written the way you would compute it.
parallel =
    auto_nomogram((R1, R2) -> 1 / (1 / R1 + 1 / R2), :R1 => (1, 10), :R2 => (1, 10), :R)
@show evaluate(parallel; R1 = 4, R2 = 6).R

# Body mass index.
bmi = auto_nomogram(
    (m, h) -> m / h^2,
    Variable(:mass; domain = (40, 150), label = "mass (kg)"),
    Variable(:height; domain = (1.4, 2.1), label = "height (m)"),
    Variable(:BMI);
    title = "Body mass index",
)
@show evaluate(bmi; mass = 70, height = 1.75).BMI

# Every chart has been checked against its formula before being returned. The relation it
# carries *is* the formula, so this measures the chart against what you wrote, not
# against anything derived along the way.
@show max_relation_residual(growth)

# Not every formula has a chart on straight scales. `separability_defect` answers the
# question directly; zero means yes.
# x + y + xy is (1+x)(1+y) - 1 in disguise, so it is separable: the defect is ~0.
println("x + y + xy:  ", separability_defect((x, y) -> x + y + x * y, (0, 5), (0, 5)))
# x² + xy + y is not.
println("x² + xy + y: ", separability_defect((x, y) -> x^2 + x * y + y, (0.5, 5), (0.5, 5)))

# And auto_nomogram says so rather than drawing something wrong.
try
    auto_nomogram((x, y) -> x^2 + x * y + y, :x => (0.5, 5), :y => (0.5, 5))
catch err
    println(sprint(showerror, err))
end

# More than two inputs: the ideal gas law. One alignment line relates only three scales,
# so this becomes two stages linked by a bare turning line -- found automatically, by
# looking for a pair of inputs the formula depends on only through one combination.
gas = auto_nomogram(
    (n, T, V) -> n * 8.314 * T / V,
    Variable(:n; domain = (0.5, 5), label = "amount (mol)"),
    Variable(:T; domain = (250, 400), label = "temperature (K)"),
    Variable(:V; domain = (0.01, 0.1), label = "volume (m³)"),
    Variable(:P; label = "pressure (Pa)");
    title = "Ideal gas",
)
@show evaluate(gas; n = 1, T = 300, V = 0.024).P

# Some formulas combine two inputs but still have no chain. In (a + b)·c, the first stage
# wants a + b on a linear scale and the second on a logarithmic one; that needs a net, and
# the error says so for each pairing it tried.
try
    auto_nomogram((a, b, c) -> (a + b) * c, :a => (1, 10), :b => (1, 10), :c => (1, 10), :d)
catch err
    println(sprint(showerror, err))
end

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)
save_svg(
    joinpath(outdir, "from_formula.svg"),
    growth;
    isopleths = [isopleth(growth; rate = 5, years = 10, label = "5% for 10 years")],
)
