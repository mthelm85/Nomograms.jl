using Pkg;
Pkg.activate(joinpath(@__DIR__, "env"); io = devnull)
using Nomograms, Roots, Printf, Random
quad = determinant_nomogram(
    Variable(:p; domain = (-12, 0)) => p -> (0.0, p),
    Variable(:q; domain = (0, 36)) => q -> (1.0, q),
    Variable(:z; domain = (0.2, 6)) => z -> (1 / (1 + z), -z^2 / (1 + z));
    relation = (p, q, z) -> z^2 + p * z + q,
)
s = quad[:z]
a, b = 1.875681, 1.888816
P = Nomograms.point(quad[:p], -(a + b));
Q = Nomograms.point(quad[:q], a * b)
g(z) = Nomograms.align_residual(s, z, P, Q)
zs = range(s.domain...; length = 4001)
r = g.(zs)
@printf(
    "residual: min %.3e  max %.3e; between the roots it dips to %.3e\n",
    minimum(r),
    maximum(r),
    g((a + b) / 2)
)
@printf(
    "sign changes on a 4001-point scan: %d\n",
    count(i -> sign(r[i]) != sign(r[i+1]), 1:(length(r)-1))
)
for n in (12, 50, 100, 200, 500)
    roots = Roots.find_zeros(g, s.domain...; no_pts = n)
    println("find_zeros no_pts = ", rpad(n, 4), " -> ", length(roots), " roots")
end
# Over the random scan: misses and cost as no_pts grows.
Random.seed!(2)
cases = Tuple{Float64,Float64}[]
for _ = 1:20000
    r1 = 0.2 + 5.8 * rand()
    r2 = 0.2 + 5.8 * rand()
    p, q = -(r1 + r2), r1 * r2
    (-12 <= p <= 0 && 0 <= q <= 36 && abs(r1 - r2) > 1e-12) && push!(cases, (r1, r2))
end
for n in (12, 50, 100, 200)
    miss = 0
    t = @elapsed for (r1, r2) in cases
        P = Nomograms.point(quad[:p], -(r1 + r2))
        Q = Nomograms.point(quad[:q], r1 * r2)
        h(z) = Nomograms.align_residual(s, z, P, Q)
        length(Roots.find_zeros(h, s.domain...; no_pts = n)) == 2 || (miss += 1)
    end
    @printf(
        "no_pts = %-4d  missed %d of %d   %.2e s per solve\n",
        n,
        miss,
        length(cases),
        t / length(cases)
    )
end
