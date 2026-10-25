# Experiment 6: the 0.6 x size x characters width estimate vs real font metrics.
using Pkg;
Pkg.activate(joinpath(@__DIR__, "env"); io = devnull)
using CairoMakie, Printf
const FT = CairoMakie.Makie.FreeTypeAbstraction
fonts = [
    ("Makie default (TeX Gyre Heros)", CairoMakie.Makie.to_font("TeX Gyre Heros Makie")),
    ("DejaVu Sans", FT.findfont("DejaVu Sans")),
    ("Arial", FT.findfont("Arial")),
]
size = 10.0
# Width of a string as laid out by FreeType: sum of glyph advances (kerning ignored).
function width(font, s, size)
    isnothing(font) && return NaN
    sum(FT.hadvance(FT.get_extent(font, c)) for c in s) * size
end
strings = [
    "1",
    "10",
    "100",
    "1000",
    "0.5",
    "1/2",
    "3/16",
    "1 3/8",
    "2 15/16",
    "-0.25",
    "R₁ (Ω)",
    "4 ∥ 6 = 2.4",
    "temperature (K)",
    "Daily Amp-Hours (Ah)",
    "Max Deflection Threshold (c)",
    "WWWW",
    "iiii",
    "1e4",
    "kg m^-2",
    "ℓ",
    "μs",
]
for (fname, font) in fonts
    isnothing(font) && (println(fname, ": not found"); continue)
    println("== ", fname, " ==")
    ratios = Float64[]
    for s in strings
        real = width(font, s, size)
        est = 0.6 * size * length(s)
        push!(ratios, real / est)
        @printf(
            "  %-30s real %6.1f  est %6.1f  real/est %.2f\n",
            repr(s),
            real,
            est,
            real / est
        )
    end
    @printf("  real/est range %.2f .. %.2f\n", minimum(ratios), maximum(ratios))
end
