using Pkg;
Pkg.activate(joinpath(@__DIR__, "env"); io = devnull)
using CairoMakie, Printf
const FT = CairoMakie.Makie.FreeTypeAbstraction
font = FT.findfont("Arial")
w(s) = sum(FT.hadvance(FT.get_extent(font, c)) for c in s) * 10
for s in [
    "MAXIMUM WATTAGE (W)",
    "HEIGHT (M)",
    "MOMENT OF INERTIA",
    "WAVELENGTH",
    "Mass (kg)",
    "Wave Height (m)",
]
    @printf("%-22s real/est %.2f\n", repr(s), w(s) / (0.6 * 10 * length(s)))
end
# a per-character Helvetica/Arial table built from the font, as the package could embed
ascii = Dict(c => FT.hadvance(FT.get_extent(font, c)) for c = Char(32):Char(126))
est2(s) = sum(get(ascii, c, 0.6) for c in s) * 10
worst = maximum(
    abs(est2(s) / w(s) - 1) for
    s in ["MAXIMUM WATTAGE (W)", "1 3/8", "temperature (K)", "Daily Amp-Hours (Ah)"]
)
@printf("per-character table: worst error %.1f%% on these strings\n", 100 * worst)
