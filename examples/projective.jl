# Reshaping a chart without changing what it computes.
#
# A nomogram asserts exactly one thing: that a point taken from each scale is collinear.
# A projective map sends straight lines to straight lines, so it preserves that and
# nothing else is sacred. The picture is free to move; the answers cannot.
#
# That turns reshaping from a risk into a design tool. If a chart is cramped at one end,
# pull it open — the arithmetic is untouched, and you do not have to re-derive anything.

using Nomograms

nom = additive_nomogram(
    Variable(:R1, ReciprocalScale(); domain = (1, 10), label = "R₁ (Ω)"),
    Variable(:R2, ReciprocalScale(); domain = (1, 10), label = "R₂ (Ω)"),
    Variable(:R, ReciprocalScale(); label = "R (Ω)");
    title = "Resistors in parallel, reshaped",
)

# Four point correspondences fix a projective map, which is the practical way to describe
# one: say where you want the corners to go. Here the bottom edge is pulled in, so the
# outer axes converge downwards.
corners = [(0, 0), (70, 0), (70, 100), (0, 100)]
skew = project(nom, corners => [(14, 0), (56, 0), (70, 100), (0, 100)])

@show evaluate(nom; R1 = 4, R2 = 6).R
@show evaluate(skew; R1 = 4, R2 = 6).R

# Not approximately the same chart — the same chart. The worst disagreement with the
# algebra across a grid of inputs is at the level of floating-point noise.
@show max_relation_residual(skew)

# A homography can also be given directly as a 3×3 matrix. The bottom row is what makes it
# projective rather than affine: with h₃₁ = h₃₂ = 0 it degenerates to an affine map, which
# is the familiar translate/scale/shear and equally safe.
@show projective_transform(corners, [(14, 0), (56, 0), (70, 100), (0, 100)])[3, :]

# The result is re-fitted to the canvas, since a projective map takes the bounding box
# somewhere arbitrary. Composing that affine fit with the projective map leaves it
# projective, so the guarantee survives.

isopleths = [isopleth(skew; R1 = 4, R2 = 6, label = "4 ∥ 6 = 2.4")]

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)
save_svg(joinpath(outdir, "projective.svg"), skew; isopleths)
