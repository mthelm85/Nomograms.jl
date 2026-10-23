@testmodule NChart begin
    using Nomograms

    speed() = ratio_nomogram(
        Variable(:distance; domain = (1, 10)),
        Variable(:time; domain = (1, 10)),
        Variable(:speed),
    )
end

@testitem "the diagonal reads the ratio" tags = [:unit] setup = [NChart] begin
    nom = NChart.speed()
    @test all(nom[:speed].domain .≈ (0.1, 10.0))
    for d in range(1, 10; length = 7), t in range(1, 10; length = 7)
        @test evaluate(nom; distance = d, time = t).speed ≈ d / t rtol = 1e-9
    end
    # Any variable may be the unknown, as always.
    @test evaluate(nom; distance = 9, speed = 2.25).time ≈ 4 rtol = 1e-9
    @test evaluate(nom; time = 4, speed = 2.25).distance ≈ 9 rtol = 1e-9
    @test relation_residual(nom, evaluate(nom; distance = 6, time = 3)) ≈ 0 atol = 1e-12
end

@testitem "the N-chart geometry is the one claimed" tags = [:unit] setup = [NChart] begin
    nom = NChart.speed()
    W, H = nom.width, nom.height
    # Outer scales run in opposite directions, anchoring the diagonal at two corners:
    # f1 upwards from the bottom-left, f2 downwards from the top-right.
    @test point(nom[:distance], 10)[2] ≈ H
    @test point(nom[:distance], 1)[2] ≈ H / 10
    @test point(nom[:time], 10)[2] ≈ 0 atol = 1e-9
    @test point(nom[:time], 1)[2] ≈ H - H / 10

    # Equal ranges give k = 1, so the ratio r sits at s = r/(r+1) along the diagonal.
    for r in (0.25, 1.0, 2.0, 5.0)
        s = r / (r + 1)
        @test all(point(nom[:speed], r) .≈ (s * W, s * H))
    end
    # The diagonal really is straight, and stops short of the corners because those are
    # the ratios 0 and infinity.
    @test length(Nomograms._path(nom[:speed])) == 2
    @test point(nom[:speed], 0.1)[1] > 0
    @test point(nom[:speed], 10)[1] < W
end

@testitem "the two outer scale factors need not match" tags = [:unit] begin
    # The textbook derivation makes them equal; here the diagonal's own scale absorbs
    # the ratio k = a1/a2, so each outer axis is fitted to the canvas independently.
    nom = ratio_nomogram(
        Variable(:a; domain = (1, 4)),
        Variable(:b; domain = (1, 20)),
        Variable(:c),
    )
    @test all(nom[:c].domain .≈ (1 / 20, 4.0))
    @test point(nom[:a], 4)[2] ≈ nom.height
    @test point(nom[:b], 20)[2] ≈ 0 atol = 1e-9
    for (x, y) in ((3, 6), (2, 15), (4, 1))
        @test evaluate(nom; a = x, b = y).c ≈ x / y rtol = 1e-9
    end
end

@testitem "transforms carry through the ratio" tags = [:unit] begin
    # f1(u)/f2(v) = f3(w) with f1 = u^2: 9/3 = 3.
    nom = ratio_nomogram(
        Variable(:u, PowerScale(2.0); domain = (1, 10)),
        Variable(:v; domain = (1, 10)),
        Variable(:w),
    )
    @test evaluate(nom; u = 3, v = 3).w ≈ 3 rtol = 1e-9
    @test evaluate(nom; u = 6, v = 4).w ≈ 9 rtol = 1e-9
end

@testitem "a ratio chart needs strictly positive scales" tags = [:unit, :validation] begin
    err = try
        ratio_nomogram(
            Variable(:a; domain = (0, 10)),
            Variable(:b; domain = (1, 10)),
            Variable(:c),
        )
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("strictly positive", err.msg)

    @test_throws ArgumentError ratio_nomogram(
        Variable(:a; domain = (1, 10)),
        Variable(:b; domain = (-5, 10)),
        Variable(:c),
    )
    @test_throws ArgumentError ratio_nomogram(
        Variable(:a),
        Variable(:b; domain = (1, 10)),
        Variable(:c),
    )
    @test_throws ArgumentError ratio_nomogram(
        Variable(:a; domain = (1, 10)),
        Variable(:b; domain = (1, 10)),
        Variable(:c);
        height = 0,
    )
end

@testitem "an inconsistent ratio domain is rejected" tags = [:unit, :validation] begin
    err = try
        ratio_nomogram(
            Variable(:a; domain = (1, 10)),
            Variable(:b; domain = (1, 10)),
            Variable(:c; domain = (1, 10)),
        )
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("inconsistent", err.msg)
    @test occursin("(0.1, 10)", err.msg)
end

@testitem "a ratio axis spanning decades gets logarithmic ticks" tags = [:unit] setup =
    [NChart] begin
    nom = NChart.speed()
    labelled = [
        t.value for t in drawing(nom).scales[3].ticks if t.level == 1 && !isnothing(t.label)
    ]
    # Linear candidates over 0.1 to 10 offer nothing between 0.1 and 2, which is most of
    # the diagonal. Decade-based ones populate the whole axis.
    @test count(<(1), labelled) >= 4
    @test 0.5 in labelled
    @test 1.0 in labelled
    # And they are exact, so they print as 0.3 rather than 0.30000000000000004.
    @test all(v -> v == parse(Float64, Nomograms.default_tickformat(v)), labelled)

    # A narrow ratio range stays linear, and the caller can always override.
    narrow = ratio_nomogram(
        Variable(:a; domain = (9, 10)),
        Variable(:b; domain = (9, 10)),
        Variable(:c),
    )
    @test Nomograms.resolved_tickstyle(narrow[:c].variable) === :linear
    forced = ratio_nomogram(
        Variable(:a; domain = (1, 10)),
        Variable(:b; domain = (1, 10)),
        Variable(:c; tickstyle = :linear),
    )
    @test Nomograms.resolved_tickstyle(forced[:c].variable) === :linear
    @test_throws ArgumentError Variable(:a; tickstyle = :loglog)
end

@testitem "an N-chart is never ambiguous" tags = [:unit] setup = [NChart] begin
    # The diagonal is straight and s = kr/(kr+1) is monotone, so there is one crossing.
    nom = NChart.speed()
    for d in (2.0, 5.0, 9.0), t in (1.5, 4.0, 8.0)
        @test length(evaluate_all(nom; distance = d, time = t)) == 1
    end
end

@testitem "the N-chart renders in both backends" tags = [:integration] setup = [NChart] begin
    using CairoMakie

    nom = NChart.speed()
    isos = [isopleth(nom; distance = 6, time = 3, label = "6 / 3 = 2")]
    out = svg(nom; isopleths = isos)
    @test occursin("id=\"scale-speed\"", out)
    @test !occursin("NaN", out)
    @test nomogramplot(nom; isopleths = isos) isa Figure
end
