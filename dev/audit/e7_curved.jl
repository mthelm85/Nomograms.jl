# Experiment 7: prototype test for the one-curved-scale (Cauchy) form.
# For F(u, v, w) = 0 with w on the curve: exists iff L = ln|F_u/F_v| = phi(w) + a(u) - b(v),
# i.e. along each level set w = c:  L_uv = A(c) w_u w_v + B(c) w_uv  with A = dB/dc.
using Pkg;
Pkg.activate(joinpath(@__DIR__, "env"); io = devnull)
using ForwardDiff, Roots, LinearAlgebra, Printf
const D = ForwardDiff

# w(u, v) solving F(u, v, w) = 0 by Newton from a float start; generic, so Duals carry exact
# derivatives of the implicit function.
function wsolve(F, u, v, w0)
    w = w0 + zero(u) + zero(v)
    for _ = 1:60
        f = F(u, v, w)
        fw = D.derivative(t -> F(u, v, t), w)
        w -= f / fw
    end
    return w
end

function cauchy_test(F, ud, vd, wd; nlev = 21, npts = 15)
    # float solution, used as the Newton start everywhere
    wf(u, v) = find_zero(t -> F(u, v, t), wd)
    W(x) = wsolve(
        F,
        x[1],
        x[2],
        wf(D.value(D.value(x[1])) |> float, D.value(D.value(x[2])) |> float),
    )
    # L as a function of (u, v) on the surface
    Lfun(x) = begin
        w = W(x)
        g = D.gradient(t -> F(t[1], t[2], t[3]), [x[1], x[2], w])
        log(abs(g[1] / g[2]))
    end
    levels =
        range(wd[1] + 0.15 * (wd[2] - wd[1]), wd[2] - 0.15 * (wd[2] - wd[1]); length = nlev)
    As, Bs, resid = Float64[], Float64[], 0.0
    for c in levels
        rows = Vector{Float64}[]
        rhs = Float64[]
        for u in range(ud[1], ud[2]; length = npts)
            # v on this level set
            local v
            try
                v = find_zero(t -> F(u, t, c), vd)
            catch
                continue
            end
            x = [u, v]
            Lxy = D.hessian(Lfun, x)[1, 2]
            gw = D.gradient(W, x)
            Hw = D.hessian(W, x)
            push!(rows, [gw[1] * gw[2], Hw[1, 2]])
            push!(rhs, Lxy)
        end
        length(rhs) < 4 && continue
        M = reduce(hcat, rows)'
        ab = M \ rhs
        r = norm(M * ab - rhs) / max(norm(rhs), 1e-12)
        norm(rhs) < 1e-10 && (r = 0.0)           # L_uv vanishes: trivially fine
        resid = max(resid, r)
        push!(As, ab[1])
        push!(Bs, ab[2])
    end
    # A = dB/dc across levels (central differences)
    h = step(levels)
    cons =
        length(Bs) >= 3 ?
        maximum(
            abs(As[i] - (Bs[i+1] - Bs[i-1]) / (2h)) /
            max(1e-9, maximum(abs, As) + maximum(abs, Bs)) for i = 2:(length(Bs)-1)
        ) : NaN
    return resid, cons
end

cases = [
    # name, F(x, y, z), domains, assignments to try: (curved var index)
    ("z = x*y (separable)", (x, y, z) -> x * y - z, (1.0, 5.0), (1.0, 5.0), (1.0, 25.0)),
    (
        "z^2 + p z + q = 0",
        (p, q, z) -> z^2 + p * z + q,
        (-7.0, -5.0),
        (1.0, 6.0),
        (3.05, 6.8),
    ),
    ("x z + y = z^3", (x, y, z) -> x * z + y - z^3, (1.0, 3.0), (1.0, 3.0), (1.5, 2.5)),
    (
        "x^2 + x y + y = z",
        (x, y, z) -> x^2 + x * y + y - z,
        (0.5, 3.0),
        (0.5, 3.0),
        (1.0, 20.0),
    ),
    (
        "sin x + e^y + x y = z",
        (x, y, z) -> sin(x) + exp(y) + x * y - z,
        (0.1, 1.0),
        (0.1, 1.0),
        (1.3, 4.5),
    ),
]
perm(F, k) = k == 3 ? F : k == 1 ? ((a, b, c) -> F(c, a, b)) : ((a, b, c) -> F(a, c, b))
names = ("x", "y", "z")
for (n, F, dx, dy, dz) in cases
    doms = (dx, dy, dz)
    print(rpad(n, 24))
    for k = 1:3
        others = [i for i = 1:3 if i != k]
        Fk = perm(F, k)
        try
            r, c = cauchy_test(Fk, doms[others[1]], doms[others[2]], doms[k])
            @printf("  %s curved: fit %.1e, A=B' %.1e", names[k], r, c)
        catch e
            @printf("  %s curved: (failed: %s)", names[k], first(sprint(showerror, e), 30))
        end
    end
    println()
end
