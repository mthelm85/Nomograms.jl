# Roots of a quadratic:
#
#     z² + pz + q = 0
#
# This one needs a curve. No arrangement of straight parallel scales can express it, so it
# is built from `determinant_nomogram`, which takes the curves directly.
#
# # Deriving the curve
#
# Put p on the line x = 0 at height p, and q on x = 1 at height q. A straight line between
# them crosses x = x₃ at height p + x₃(q - p), so a third scale at (x₃(z), y₃(z)) is met
# exactly when
#
#     (1 - x₃)·p + x₃·q = y₃
#
# We want that to say the same thing as z·p + 1·q + z² = 0, so match coefficients:
#
#     (1 - x₃)/z = x₃/1 = -y₃/z²
#
# The first equality gives x₃ = 1/(1 + z); substituting into the second gives
# y₃ = -z²/(1 + z). That is the curve, and it is all the construction needs.
#
# # Two roots, two crossings
#
# A quadratic has two roots, and the nomogram shows why: there is one alignment line, and
# it meets the curved scale twice. `evaluate` will not choose between them -- it throws and
# names both -- so use `evaluate_all`.

using Nomograms

nom = determinant_nomogram(
    Variable(:p; domain = (-6, 0), label = "p") => p -> (0.0, p),
    Variable(:q; domain = (0, 6), label = "q") => q -> (1.0, q),
    Variable(:z; domain = (0.2, 5), label = "z") => z -> (1 / (1 + z), -z^2 / (1 + z));
    relation = (p, q, z) -> z^2 + p * z + q,
    title = "z² + pz + q = 0",
)

# The curves were written in their own coordinates -- x from 0 to 1, y from -6 to 6 -- and
# fitted to the canvas. Collinearity survives any affine map, so the fit cannot change what
# the chart computes.
@show nom[:z].domain

# Does the geometry really encode the formula? This solves on the chart and checks the
# answers against the algebra, which is the question you actually want answered after
# deriving a determinant form by hand.
@show max_relation_residual(nom)

# z² - 3z + 2 = (z - 1)(z - 2)
@show [sol.z for sol in evaluate_all(nom; p = -3, q = 2)]

# z² - 5z + 4 = (z - 1)(z - 4)
@show [sol.z for sol in evaluate_all(nom; p = -5, q = 4)]

# Both isopleths are the *same* straight line -- that is the point -- so they are drawn on
# top of one another, with a marker at each crossing of the curve.
isopleths = [
    isopleth(nom, sol; label = "z = $(round(sol.z; digits = 3))") for
    sol in evaluate_all(nom; p = -3, q = 2)
]

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)
save_svg(joinpath(outdir, "quadratic.svg"), nom; isopleths)
