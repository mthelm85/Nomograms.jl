@testmodule Comp begin
    using Nomograms

    # V = L * W * H, split as t = L*W then V = t*H. Both stages are multiplicative, so
    # both place the shared `t` logarithmically and the two scales can coincide.
    function volume(; blank = false, width2 = 130)
        tk = blank ? Float64[] : nothing
        s1 = multiplicative_nomogram(
            Variable(:L; domain = (1, 10)),
            Variable(:W; domain = (1, 10)),
            Variable(:t; label = "", ticks = tk),
        )
        s2 = multiplicative_nomogram(
            Variable(:t; domain = s1[:t].domain, label = "", ticks = tk),
            Variable(:H; domain = (1, 10)),
            Variable(:V);
            width = width2,
        )
        return (s1, s2, compound(s1, s2))
    end
end

@testitem "a compound chart solves in any direction" tags = [:unit] setup = [Comp] begin
    _, _, nom = Comp.volume()
    @test Set(variablenames(nom)) == Set([:L, :W, :t, :H, :V])

    for (L, W, H) in ((3, 4, 2), (2, 2, 5), (9, 1, 1), (10, 10, 10), (1.5, 7.5, 3.25))
        got = evaluate(nom; L, W, H)
        @test got.t ≈ L * W rtol = 1e-8
        @test got.V ≈ L * W * H rtol = 1e-8
    end
    # Backwards through the second stage...
    @test evaluate(nom; L = 3, W = 4, V = 24).H ≈ 2 rtol = 1e-8
    # ...and backwards through the first, with the intermediate supplied instead.
    @test evaluate(nom; L = 3, t = 12, H = 2).W ≈ 4 rtol = 1e-8
    # But knowing the intermediate says nothing about how it was reached: t = 12 leaves
    # L and W with only their product fixed, so this is genuinely under-determined.
    @test_throws ArgumentError evaluate(nom; t = 12, H = 2, V = 24)
end

@testitem "the shared scale is drawn once, and is really shared" tags = [:unit] setup =
    [Comp] begin
    s1, s2, nom = Comp.volume()
    # Five scales from two stages of three: `t` appears once.
    @test length(nom.scales) == 5
    @test count(n -> n === :t, variablenames(nom)) == 1

    # Both stages' own copies of :t now sit exactly on the drawn one.
    for stage in nom.stages, v in range(1, 100; length = 9)
        @test all(point(stage[:t], v) .≈ point(nom[:t], v))
    end
    # Left to right, so the scale list matches the picture.
    xs = [point(s, (s.domain[1] + s.domain[2]) / 2)[1] for s in nom.scales]
    @test issorted(xs)
end

@testitem "alignment cannot change what a stage computes" tags = [:unit] setup = [Comp] begin
    # A similarity is affine, so each stage must still answer exactly as it did alone.
    s1, s2, nom = Comp.volume()
    for L in range(1, 10; length = 7), W in range(1, 10; length = 7)
        @test evaluate(nom.stages[1]; L, W).t ≈ evaluate(s1; L, W).t rtol = 1e-9
    end
    @test max_relation_residual(nom.stages[1]) < 1e-9
    @test max_relation_residual(nom.stages[2]) < 1e-9
end

@testitem "stages that cannot share a scale are rejected" tags = [:unit, :validation] begin
    # The first stage puts :t on a linear scale, the second on a logarithmic one. No
    # amount of moving the stages reconciles those, and a chart drawn anyway would lie.
    s1 = additive_nomogram(
        Variable(:a; domain = (1, 10)),
        Variable(:b; domain = (1, 10)),
        Variable(:t),
    )
    s2 = multiplicative_nomogram(
        Variable(:t; domain = s1[:t].domain),
        Variable(:c; domain = (1, 10)),
        Variable(:d),
    )
    err = try
        compound(s1, s2)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("on top of one another", err.msg)

    # A domain mismatch is reported separately, with the domain to use.
    s3 = multiplicative_nomogram(
        Variable(:t; domain = (1, 50)),
        Variable(:c; domain = (1, 10)),
        Variable(:d),
    )
    s0 = multiplicative_nomogram(
        Variable(:a; domain = (1, 10)),
        Variable(:b; domain = (1, 10)),
        Variable(:t),
    )
    err = try
        compound(s0, s3)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("disagree about :t", err.msg)
    @test occursin("(1, 100)", err.msg)
end

@testitem "compound argument errors" tags = [:unit, :validation] setup = [Comp] begin
    s1, s2, nom = Comp.volume()
    @test_throws ArgumentError compound(s1)

    # Nothing in common to pivot about.
    lone = multiplicative_nomogram(
        Variable(:x; domain = (1, 10)),
        Variable(:y; domain = (1, 10)),
        Variable(:z),
    )
    err = try
        compound(s1, lone)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("shares no variable", err.msg)

    @test_throws ArgumentError evaluate(nom; L = 1, W = 1, t = 1, H = 1, V = 1)
    @test_throws ArgumentError evaluate(nom; L = 3, q = 1)

    # Not enough to close the chain: the first stage still has two unknowns.
    err = try
        evaluate(nom; L = 3, H = 2)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("stayed unknown", err.msg)
end

@testitem "a reading is one group across both stages" tags = [:unit] setup = [Comp] begin
    _, _, nom = Comp.volume()
    isos = isopleth(nom; L = 3, W = 4, H = 2, label = "3x4x2")
    @test length(isos) == 2
    @test allequal(i.group for i in isos)
    # The label goes where the answer is read.
    @test isos[1].label == ""
    @test isos[2].label == "3x4x2"
    @test isos[1][:t] ≈ 12 rtol = 1e-8
    @test isos[2][:V] ≈ 24 rtol = 1e-8

    d = drawing(nom; isopleths = isos)
    @test length(d.isopleths) == 2
    @test all(id -> id.group == isos[1].group, d.isopleths)

    # From a solution you already have.
    sol = evaluate(nom; L = 3, W = 4, H = 2)
    @test length(isopleth(nom, sol; group = 3)) == 2
end

@testitem "grouped isopleths share a colour" tags = [:unit] begin
    colors = Nomograms.DEFAULT_ISOPLETH_COLORS
    dashes = Nomograms.DEFAULT_ISOPLETH_DASHES
    # Two segments of one reading must not be painted as two different series.
    @test Nomograms.isopleth_style(2, colors, dashes) ==
          Nomograms.isopleth_style(2, colors, dashes)
    @test Nomograms.isopleth_style(1, colors, dashes) !=
          Nomograms.isopleth_style(2, colors, dashes)
end

@testitem "an empty tick list draws a bare turning line" tags = [:unit] setup = [Comp] begin
    _, _, plain = Comp.volume(; blank = false)
    _, _, blank = Comp.volume(; blank = true)
    tick_count(nom) = length(only(sd for sd in drawing(nom).scales if sd.name === :t).ticks)
    @test tick_count(plain) > 0
    # Endpoint ticks would otherwise survive an explicit "no ticks".
    @test tick_count(blank) == 0
end

@testitem "compound charts render in both backends" tags = [:integration] setup = [Comp] begin
    using CairoMakie

    _, _, nom = Comp.volume(; blank = true)
    isos = isopleth(nom; L = 3, W = 4, H = 2, label = "3x4x2 = 24")
    out = svg(nom; isopleths = isos)
    @test occursin("id=\"scale-V\"", out)
    @test occursin("id=\"isopleth-2\"", out)
    @test !occursin("NaN", out)
    @test nomogramplot(nom; isopleths = isos) isa Figure
end
