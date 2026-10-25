using Pkg;
Pkg.activate(joinpath(@__DIR__, "env"); io = devnull)
using Nomograms, Printf, Random
quad = determinant_nomogram(
    Variable(:p; domain = (-12, 0)) => p -> (0.0, p),
    Variable(:q; domain = (0, 36)) => q -> (1.0, q),
    Variable(:z; domain = (0.2, 6)) => z -> (1 / (1 + z), -z^2 / (1 + z));
    relation = (p, q, z) -> z^2 + p * z + q,
)
found(p, q) =
    try
        sort([s.z for s in evaluate_all(quad; p, q)])
    catch
        Float64[]
    end
Random.seed!(2)
misses = Tuple{Float64,Float64}[]
for _ = 1:20000
    r1 = 0.2 + 5.8 * rand()
    r2 = 0.2 + 5.8 * rand()
    abs(r1 - r2) < 1e-12 && continue
    p, q = -(r1 + r2), r1 * r2
    (-12 <= p <= 0 && 0 <= q <= 36) || continue
    length(found(p, q)) == 2 || push!(misses, minmax(r1, r2))
end
println(length(misses), " misses in 20000")
for (a, b) in misses[1:min(end, 8)]
    @printf(
        "roots %.6f, %.6f  gap %.2e  found %s\n",
        a,
        b,
        b - a,
        string(round.(found(-(a + b), a * b); digits = 6))
    )
end
# Look at the residual the solver searches for one miss: the signed distance from the line.
if !isempty(misses)
    a, b = misses[1]
    p, q = -(a + b), a * b
    P = Nomograms.point(quad[:p], p)
    Q = Nomograms.point(quad[:q], q)
    s = quad[:z]
    res(z) = Nomograms.align_residual(Nomograms.point(s, z), P, Q)
    zs = range(s.domain...; length = 2001)
    r = res.(zs)
    @printf(
        "residual range on the scale: %.3e .. %.3e; at the roots: %.2e, %.2e\n",
        minimum(r),
        maximum(r),
        res(a),
        res(b)
    )
    i = argmin(abs.(r))
    @printf("smallest |residual| on a 2001-point grid at z = %.4f\n", zs[i])
end
