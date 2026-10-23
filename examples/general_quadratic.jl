# The general quadratic:
#
#     a·x² + b·x + c = 0
#
# Four variables. Chaining two stages through a pivot (see `box_volume.jl`) will not do it:
# a pivot works only when two of the variables collapse into a single intermediate, and
# here they do not. This wants a **net** — a two-parameter mesh of curves occupying one of
# the three points an alignment line crosses, so that one point carries two variables.
#
# # Deriving the net
#
# Put `a` on x = 0 at height `a`, and `b` on x = 1 at height `b`. A line between them meets
# x = x₃ at height (1-x₃)·a + x₃·b, so a net point (x₃, y₃) is crossed exactly when
#
#     α·a + β·b = γ        with   x₃ = β/(α+β),  y₃ = γ/(α+β)
#
# Writing the quadratic as a·x² + b·x = -c gives α = x², β = x, γ = -c, so
#
#     x₃ = x/(x² + x) = 1/(1 + x)        y₃ = -c/(x² + x) = -c/(x·(1 + x))
#
# Both `x` and `c` move the point, and independently — which is what makes this a genuine
# net rather than a pivot in disguise.
#
# # Reading it
#
# Lay a line from `a` to `b`. Where it crosses the strand for your `c`, read off the `x`
# strand through that point. There are generally two such crossings, because a quadratic
# has two roots.

using Nomograms

nom = net_nomogram(
    Net(
        Variable(:x; domain = (0.5, 5), label = "x"),
        Variable(:c; domain = (0, 4), label = "c"),
        (x, c) -> (1 / (1 + x), -c / (x * (1 + x))),
    ),
    Variable(:a; domain = (0.5, 3), label = "a") => a -> (0.0, a),
    Variable(:b; domain = (-6, 0), label = "b") => b -> (1.0, b);
    relation = (x, c, a, b) -> a * x^2 + b * x + c,
    title = "a·x² + b·x + c = 0",
)

# Does the mesh really encode the formula? This sweeps the chart and checks each answer
# against the algebra, which is the question worth asking after deriving a net by hand.
@show max_relation_residual(nom)

# x² - 3x + 2 = (x-1)(x-2)
@show [sol.x for sol in evaluate_all(nom; a = 1, b = -3, c = 2)]

# Halving a and doubling nothing else is a different equation with the same roots.
@show [sol.x for sol in evaluate_all(nom; a = 2, b = -6, c = 4)]

# Any of the four may be the unknown. With both net variables known the net is a single
# point and this is the ordinary three-point problem.
@show evaluate(nom; x = 2, c = 2, b = -3).a
@show evaluate(nom; x = 2, a = 1, b = -3).c

isopleths = [
    isopleth(nom, sol; label = "x = $(round(sol.x; digits = 2))") for
    sol in evaluate_all(nom; a = 1, b = -3, c = 2)
]

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)
save_svg(joinpath(outdir, "general_quadratic.svg"), nom; isopleths)
