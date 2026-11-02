# Regenerates the pictures the README shows. Run from the repository root with
#
#     julia --project=. dev/readme_images.jl

using Nomograms

const OUT = joinpath(@__DIR__, "..", "docs", "src", "assets")
mkpath(OUT)

parallel = additive_nomogram(
    Variable(:R1, ReciprocalScale(); domain = (1, 10), label = "R₁ (Ω)"),
    Variable(:R2, ReciprocalScale(); domain = (1, 10), label = "R₂ (Ω)"),
    Variable(:R, ReciprocalScale(); label = "R (Ω)");
    title = "Resistors in parallel",
)
save_svg(
    joinpath(OUT, "readme-parallel.svg"),
    parallel;
    isopleths = [isopleth(parallel; R1 = 4, R2 = 6, label = "2.4 Ω")],
)

z(u, v) = (1 + complex(u, v)) / (1 - complex(u, v))
major = [0.2, 0.5, 1, 2, 5, 10, 20]
minor = [0.1, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1.2, 1.4, 1.6, 1.8, 3, 4]
smith = contour_chart(
    :u => (-1, 1),
    :v => (-1, 1),
    Variable(
        :r;
        domain = (0, 50),
        ticks = major,
        minorticks = minor,
        label = "resistance",
    ) => (u, v) -> real(z(u, v)),
    Variable(
        :x;
        domain = (-50, 50),
        ticks = [0; major; -major],
        minorticks = [minor; -minor],
        label = "reactance",
    ) => (u, v) -> imag(z(u, v));
    axes = false,
    region = (u, v) -> u^2 + v^2 <= 1,
    rings = [
        Variable(
            :angle;
            domain = (-180, 180),
            ticks = -150:30:180,
            minorticks = -180:10:180,
            label = "angle of reflection coefficient (°)",
        ) => θ -> (1.06cosd(θ), 1.06sind(θ)),
    ],
    labelat = (r = :x => 0, x = :r => 0),
    title = "A Smith chart",
)
save_svg(
    joinpath(OUT, "readme-smith.svg"),
    smith;
    crossings = [crossing(smith; r = 1, x = 1, label = "1 + j")],
)
