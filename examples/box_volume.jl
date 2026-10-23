# Four variables:
#
#     V = L × W × H
#
# One alignment line can only relate three scales, so four variables need two stages
# sharing an intermediate:
#
#     t = L × W        then        V = t × H
#
# The reader lays a line from L to W, notes where it crosses the shared scale, then pivots
# about that point and reads a line to H. Hence the usual name for the shared scale: the
# turning line.
#
# # The one constraint
#
# A shared scale has to be the *same* scale in both stages -- same values, same positions.
# Both stages here are multiplicative, so both place `t` logarithmically and they match.
# Splitting `V = (L + W) × H` the same way would not work: the first stage would put `t`
# on a linear scale and the second on a logarithmic one, and no amount of moving the
# stages can reconcile those. `compound` checks and says so rather than drawing a chart
# that lies.

using Nomograms

# `t` is given an empty tick list, which draws it as a bare turning line. Drop that and it
# carries its own scale, which is useful if the intermediate means something.
stage1 = multiplicative_nomogram(
    Variable(:L; domain = (1, 10), label = "L (m)"),
    Variable(:W; domain = (1, 10), label = "W (m)"),
    Variable(:t; label = "", ticks = Float64[]),
)

# The second stage is built with the domain the first one derived -- that is what makes
# the two `t` scales line up. Its wider canvas keeps its scales clear of the first
# stage's rather than interleaving with them.
stage2 = multiplicative_nomogram(
    Variable(:t; domain = stage1[:t].domain, label = "", ticks = Float64[]),
    Variable(:H; domain = (1, 10), label = "H (m)"),
    Variable(:V; label = "V (m³)");
    width = 130,
)

nom = compound(stage1, stage2; title = "V = L × W × H")

@show stage1[:t].domain
@show evaluate(nom; L = 3, W = 4, H = 2)

# Any direction, as long as each stage ends up with one unknown.
@show evaluate(nom; L = 3, W = 4, V = 24).H
@show evaluate(nom; L = 3, t = 12, H = 2).W

# One reading, two segments: they share a `group`, so both are drawn in the same colour.
isopleths = isopleth(nom; L = 3, W = 4, H = 2, label = "3 × 4 × 2 = 24")

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)
save_svg(joinpath(outdir, "box_volume.svg"), nom; isopleths)
