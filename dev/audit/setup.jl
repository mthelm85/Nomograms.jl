# Creates dev/audit/env with Nomograms (this checkout), ForwardDiff, CairoMakie, Unitful
# and Roots. Run once: julia dev/audit/setup.jl; then julia dev/audit/e1_thresholds.jl etc.
using Pkg
Pkg.activate(joinpath(@__DIR__, "env"))
Pkg.develop(path = joinpath(@__DIR__, "..", ".."))
Pkg.add(["ForwardDiff", "CairoMakie", "Unitful", "Roots"])
using ForwardDiff, CairoMakie, Nomograms
println("ok ", pkgversion(ForwardDiff), " ", pkgversion(CairoMakie))
