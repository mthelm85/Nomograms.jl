# Compound growth:
#
#     growth = (1 + r)^n
#
# Four-variable money problems (PV, FV, rate, years) are *not* nomographable on three
# parallel scales, because PV and FV only ever appear as their ratio. Chart the ratio and
# the relation separates:
#
#     log(growth) = n * log(1 + r)
#
# which is a product of a function of `n` and a function of `r`, so it is exactly what
# `multiplicative_nomogram` is for. Each variable keeps its own transform underneath the
# logarithm that the construction adds:
#
#   years   f(n) = n                     -- a plain linear transform
#   rate    f(r) = log(1 + r/100)        -- a custom one, with its inverse
#   growth  f(g) = log(g)                -- a natural logarithm
#
# Reading the chart: a line from a rate to a number of years crosses `growth` at the
# multiple your money becomes.

using Nomograms

rate_transform = FunctionScale(r -> log(1 + r / 100), y -> 100 * expm1(y))

nom = multiplicative_nomogram(
    Variable(:years, LinearScale(); domain = (1, 30), label = "years"),
    Variable(:rate, rate_transform; domain = (1, 20), label = "rate (% p.a.)"),
    Variable(
        :growth,
        LogScale(ℯ);
        label = "growth multiple",
        # Position on this axis goes as log(log(growth)), which compresses the top end
        # savagely: automatic decade ticks would label almost nothing above 10. Explicit
        # ticks fix that. Anything that will not fit is still dropped, so the list can be
        # generous.
        ticks = [1.02, 1.03, 1.05, 1.1, 1.25, 1.5, 2, 3, 5, 10, 20, 50, 100, 200],
    );
    title = "Compound growth: (1 + r)ⁿ",
)

@show nom[:growth].domain                          # about (1.01, 237)
@show evaluate(nom; rate = 5, years = 10).growth   # 1.05^10 = 1.6289...

# Inverting the chart answers the question people actually ask: at 7%, how long until
# the money doubles? (The rule of 72 says about 10.3 years.)
@show evaluate(nom; rate = 7, growth = 2).years

isopleths = [
    isopleth(nom; rate = 5, years = 10, label = "5% for 10y"),
    isopleth(nom; rate = 7, growth = 2, label = "7% doubles"),
]

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)
save_svg(joinpath(outdir, "compound_growth.svg"), nom; isopleths)
