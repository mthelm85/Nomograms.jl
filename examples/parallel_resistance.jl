# Two resistors in parallel:
#
#     1/R = 1/R1 + 1/R2
#
# This is already in additive form -- the transform is simply 1/x -- so it needs no
# rearranging. Note that `R` declares no domain: the range it can reach is fixed by the
# other two scales, and the construction derives it (here, 0.5 to 5 ohm).

using Nomograms

nom = additive_nomogram(
    Variable(:R1, ReciprocalScale(); domain = (1, 10), label = "R₁ (Ω)"),
    Variable(:R2, ReciprocalScale(); domain = (1, 10), label = "R₂ (Ω)"),
    Variable(:R, ReciprocalScale(); label = "R (Ω)");
    title = "Resistors in parallel",
)

@show nom[:R].domain                      # (0.5, 5.0)
@show evaluate(nom; R1 = 4, R2 = 6)       # R = 2.4
@show evaluate(nom; R1 = 4, R = 2.4)      # ...and backwards: R2 = 6

isopleths = [
    isopleth(nom; R1 = 4, R2 = 6, label = "4 ∥ 6 = 2.4"),
    isopleth(nom; R1 = 1.5, R2 = 9, label = "1.5 ∥ 9 ≈ 1.29"),
]

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)

# The SVG backend is built in and needs no plotting package.
save_svg(joinpath(outdir, "parallel_resistance.svg"), nom; isopleths)

# With a Makie backend loaded you get a Figure instead, and Makie's `save`.
#
#     using CairoMakie
#     fig = nomogramplot(nom; isopleths)
#     save(joinpath(outdir, "parallel_resistance.pdf"), fig)
