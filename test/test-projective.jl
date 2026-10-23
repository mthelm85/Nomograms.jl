@testmodule Proj begin
    using Nomograms

    chart() = additive_nomogram(
        Variable(:a; domain = (0, 10)),
        Variable(:b; domain = (0, 10)),
        Variable(:s),
    )

    # A keystone: the bottom edge pulled in, the top left alone. h31/h32 come out nonzero,
    # so this is a genuine perspective rather than an affine map in disguise.
    const SRC = [(0, 0), (70, 0), (70, 100), (0, 100)]
    const DST = [(10, 0), (60, 0), (70, 100), (0, 100)]
end

@testitem "a four-point map really sends src to dst" tags = [:unit] setup = [Proj] begin
    H = projective_transform(Proj.SRC, Proj.DST)
    for (s, d) in zip(Proj.SRC, Proj.DST)
        x, y = float(s[1]), float(s[2])
        w = H[3, 1] * x + H[3, 2] * y + H[3, 3]
        got = (
            (H[1, 1] * x + H[1, 2] * y + H[1, 3]) / w,
            (H[2, 1] * x + H[2, 2] * y + H[2, 3]) / w,
        )
        @test got[1] ≈ d[1] atol = 1e-9
        @test got[2] ≈ d[2] atol = 1e-9
    end
    # Genuinely projective, not an affine map wearing a 3x3 coat.
    @test !(abs(H[3, 1]) < 1e-12 && abs(H[3, 2]) < 1e-12)
end

@testitem "projection changes the picture, not the answers" tags = [:unit] setup = [Proj] begin
    nom = Proj.chart()
    skew = project(nom, Proj.SRC => Proj.DST)

    # This is the whole guarantee: a projective map sends lines to lines, so collinearity
    # -- the only thing a nomogram asserts -- survives it.
    errs = [
        abs(evaluate(skew; a, b).s - (a + b)) for
        a in range(0, 10; length = 11), b in range(0, 10; length = 11)
    ]
    @test maximum(errs) < 1e-9
    @test max_relation_residual(skew) < 1e-9

    # But the geometry did move.
    @test point(skew[:a], 0) != point(nom[:a], 0)
    @test skew.relation === nom.relation
end

@testitem "the identity is a no-op" tags = [:unit] setup = [Proj] begin
    nom = Proj.chart()
    same = project(nom, [1.0 0 0; 0 1.0 0; 0 0 1.0])
    for v in (0.0, 3.3, 10.0)
        @test all(point(same[:a], v) .≈ point(nom[:a], v))
        @test all(point(same[:s], v) .≈ point(nom[:s], v))
    end
end

@testitem "a projected chart still fills its canvas" tags = [:unit] setup = [Proj] begin
    skew = project(Proj.chart(), Proj.SRC => Proj.DST; width = 120, height = 60)
    @test skew.width == 120 && skew.height == 60
    xs = Float64[]
    ys = Float64[]
    for s in skew.scales, t in range(s.domain...; length = 33)
        p = point(s, t)
        push!(xs, p[1])
        push!(ys, p[2])
    end
    @test minimum(xs) ≈ 0 atol = 1e-6
    @test maximum(xs) ≈ 120 atol = 1e-6
    @test minimum(ys) ≈ 0 atol = 1e-6
    @test maximum(ys) ≈ 60 atol = 1e-6
end

@testitem "projection works on a derived result domain" tags = [:unit] begin
    # :s leaves its domain to be derived, so the Variable has none and only its Scale
    # knows the range -- which is what the fit has to read.
    nom = additive_nomogram(
        Variable(:a; domain = (0, 10)),
        Variable(:b; domain = (0, 20)),
        Variable(:s),
    )
    @test isnothing(nom[:s].variable.domain)
    skew = project(
        nom,
        [(0, 0), (70, 0), (70, 100), (0, 100)] =>
            [(5, 0), (65, 0), (70, 100), (0, 100)],
    )
    @test all(skew[:s].domain .≈ (0.0, 30.0))
    @test evaluate(skew; a = 3, b = 4).s ≈ 7 atol = 1e-9
end

@testitem "projection composes and survives a curved chart" tags = [:unit] begin
    quad = determinant_nomogram(
        Variable(:p; domain = (-6, 0)) => p -> (0.0, p),
        Variable(:q; domain = (0, 6)) => q -> (1.0, q),
        Variable(:z; domain = (0.2, 5)) => z -> (1 / (1 + z), -z^2 / (1 + z));
        relation = (p, q, z) -> z^2 + p * z + q,
    )
    src = [(0, 0), (70, 0), (70, 100), (0, 100)]
    once = project(quad, src => [(8, 0), (62, 0), (70, 100), (0, 100)])
    twice = project(once, src => [(0, 0), (70, 0), (62, 100), (8, 100)])
    for nom in (once, twice)
        @test max_relation_residual(nom) < 1e-7
        got = [s.z for s in evaluate_all(nom; p = -3, q = 2)]
        @test length(got) == 2
        @test got[1] ≈ 1 atol = 1e-6
        @test got[2] ≈ 2 atol = 1e-6
    end
end

@testitem "projective argument errors" tags = [:unit, :validation] setup = [Proj] begin
    nom = Proj.chart()
    @test_throws ArgumentError project(nom, [1.0 0; 0 1.0])
    @test_throws ArgumentError project(nom, [1.0 1 1; 2 2 2; 3 3 3])   # singular
    @test_throws ArgumentError project(nom, [1.0 0 0; 0 1.0 0; 0 0 1.0]; width = 0)
    @test_throws ArgumentError project(nom, [1.0 0 0; 0 1.0 0; 0 NaN 1.0])

    # The vanishing line w = 1 - y/50 = 0 runs straight through the chart.
    err = try
        project(nom, [1.0 0 0; 0 1.0 0; 0 -1/50 1.0])
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("infinity", err.msg)

    # Four points with three collinear do not fix a map.
    @test_throws ArgumentError projective_transform(
        [(0, 0), (1, 1), (2, 2), (5, 0)],
        Proj.DST,
    )
    @test_throws ArgumentError projective_transform([(0, 0), (1, 1)], [(0, 0), (1, 1)])
end

@testitem "a reshaped chart renders" tags = [:integration] setup = [Proj] begin
    using CairoMakie

    skew = project(Proj.chart(), Proj.SRC => Proj.DST)
    isos = [isopleth(skew; a = 3, b = 4, label = "3 + 4")]
    out = svg(skew; isopleths = isos)
    @test !occursin("NaN", out)
    @test nomogramplot(skew; isopleths = isos) isa Figure
end
