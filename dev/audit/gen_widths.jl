using Pkg;
Pkg.activate(joinpath(@__DIR__, "env"); io = devnull)
using CairoMakie
const FT = CairoMakie.Makie.FreeTypeAbstraction
font = FT.findfont("Arial")
extra = collect("°±×÷µ²³¹½¼¾−–—…·•→←↑↓↔∥∠≈≠≤≥∞√∑∫πΩΔαβγδεθλμσφωℓ′″€£₀₁₂₃₄₅₆₇₈₉⁰⁴⁵⁶⁷⁸⁹⁺⁻")
chars = [Char(c) for c = 32:126]
append!(chars, extra)
io = IOBuffer()
println(
    io,
    "# Advance widths in em of Helvetica's metrics (as Arial, which matches them), for the",
)
println(io, "# characters this package and its units print. Generated from the font by")
println(
    io,
    "# dev/audit/gen_widths.jl; anything else falls back to 0.6 em, or a full em when wide.",
)
println(io, "const _ADVANCE = Dict{Char,Float64}(")
line = "   "
for c in chars
    w = FT.hadvance(FT.get_extent(font, c))
    item = " $(repr(c)) => $(round(w; digits = 4)),"
    if length(line) + length(item) > 90
        println(io, line)
        global line = "   "
    end
    global line *= item
end
println(io, line)
println(io, ")")
write(joinpath(@__DIR__, "widths.jl"), String(take!(io)))
println("wrote ", length(chars), " widths")
