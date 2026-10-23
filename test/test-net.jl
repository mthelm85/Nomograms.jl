@testmodule NetFix begin
    using Nomograms

    # a*x^2 + b*x + c = 0. Written a*x^2 + b*x = -c the coefficients are alpha = x^2,
    # beta = x, gamma = -c, so x3 = 1/(1+x) and y3 = -c/(x*(1+x)).
    quadratic(; kwargs...) = net_nomogram(
        Net(
            Variable(:x; domain = (0.5, 5)),
            Variable(:c; domain = (0, 4)),
            (x, c) -> (1 / (1 + x), -c / (x * (1 + x))),
        ),
        Variable(:a; domain = (0.5, 3)) => a -> (0.0, a),
        Variable(:b; domain = (-6, 0)) => b -> (1.0, b);
        relation = (x, c, a, b) -> a * x^2 + b * x + c,
        kwargs...,
    )
end

@testitem "a net chart solves whichever variable is unknown" tags = [:unit] setup = [NetFix] begin
    nom = NetFix.quadratic()
    @test variablenames(nom) == [:x, :c, :a, :b]

    # Both net variables unknown is not offered; one at a time is.
    for (a, b, c, roots) in
        ((1, -3, 2, (1.0, 2.0)), (2, -6, 4, (1.0, 2.0)), (0.5, -3, 4, (2.0, 4.0)))
        got = [s.x for s in evaluate_all(nom; a, b, c)]
        @test length(got) == length(roots)
        for (g, r) in zip(got, roots)
            @test g ≈ r atol = 1e-6
        end
    end
    # With both net variables known the net is a single point: the ordinary problem.
    @test evaluate(nom; x = 2, c = 2, b = -3).a ≈ 1 atol = 1e-6
    @test evaluate(nom; x = 2, c = 2, a = 1).b ≈ -3 atol = 1e-6
    # And solving for the other net variable.
    @test evaluate(nom; x = 2, a = 1, b = -3).c ≈ 2 atol = 1e-6

    @test max_relation_residual(nom) < 1e-8
end

@testitem "two roots are two crossings of one strand" tags = [:unit] setup = [NetFix] begin
    nom = NetFix.quadratic()
    # The line is fixed by a and b; the c strand is a curve, so it is crossed twice.
    @test length(evaluate_all(nom; a = 1, b = -3, c = 2)) == 2
    err = try
        evaluate(nom; a = 1, b = -3, c = 2)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("crosses the net at 2 values of :x", err.msg)
    @test occursin("evaluate_all", err.msg)

    # A line that misses the strand entirely is reported, not fudged.
    @test_throws ArgumentError evaluate_all(nom; a = 3, b = -0.5, c = 4)
end

@testitem "the net is drawn as a light mesh" tags = [:unit] setup = [NetFix] begin
    nom = NetFix.quadratic()
    grid = [s for s in nom.scales if s.kind === :grid]
    plain = [s for s in nom.scales if s.kind === :scale]
    @test length(plain) == 2
    @test length(grid) > 4
    # Strands carry no ticks of their own; they are labelled once at an end.
    d = drawing(nom)
    for sd in d.scales
        sd.kind === :grid || continue
        @test isempty(sd.ticks)
        @test !isempty(sd.label)
    end
    @test any(sd -> sd.kind === :scale && !isempty(sd.ticks), d.scales)

    # An ordinary scale defaults to :scale, so nothing else changes.
    v = Variable(:q; domain = (0, 1))
    @test Scale(v, (0.0, 1.0), t -> (0.0, t), :left).kind === :scale
end

@testitem "net variables have no scale of their own" tags = [:unit, :validation] setup =
    [NetFix] begin
    nom = NetFix.quadratic()
    @test nom[:a] isa Scale
    err = try
        nom[:x]
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("carried by this chart's net", err.msg)
    @test_throws KeyError nom[:nope]
end

@testitem "net_nomogram validates its arguments" tags = [:unit, :validation] begin
    good = Variable(:a; domain = (0, 1)) => a -> (0.0, a)
    net =
        Net(Variable(:u; domain = (0, 1)), Variable(:v; domain = (0, 1)), (u, v) -> (u, v))
    @test_throws ArgumentError net_nomogram(
        Net(Variable(:u), Variable(:v; domain = (0, 1)), (u, v) -> (u, v)),
        good,
        Variable(:b; domain = (0, 1)) => b -> (1.0, b),
    )
    @test_throws ArgumentError net_nomogram(net, good, Variable(:b) => b -> (1.0, b))
    @test_throws ArgumentError net_nomogram(net, :nope => a -> (0.0, a), good)
end

@testitem "too few or too many values" tags = [:unit, :validation] setup = [NetFix] begin
    nom = NetFix.quadratic()
    @test_throws ArgumentError evaluate(nom; a = 1, b = -3)
    @test_throws ArgumentError evaluate(nom; a = 1, b = -3, c = 2, x = 1)
    @test_throws ArgumentError evaluate(nom; a = 1, b = -3, zzz = 2)
end

@testitem "titleend picks which end a title sits off" tags = [:unit] begin
    # A diagonal scale, so the two ends are at different heights and :auto has a view.
    mk(te) = Scale(
        Variable(:d; domain = (0, 100), label = "d", titleend = te),
        (0.0, 100.0),
        t -> (t, t),
        :right,
    )
    flank(x, n) = Scale(Variable(n; domain = (0, 100)), (0.0, 100.0), y -> (x, y), :left)
    pos(te) = begin
        nom = Nomogram([flank(-50.0, :u), flank(150.0, :v), mk(te)]; width = 100, height = 100)
        drawing(nom).scales[3].titlepos
    end
    @test pos(:auto)[2] > 50          # the higher end
    @test pos(:stop)[2] > 50
    @test pos(:start)[2] < 50         # forced to the low end
    @test_throws ArgumentError Variable(:a; titleend = :middle)
end

@testitem "net charts render in both backends" tags = [:integration] setup = [NetFix] begin
    using CairoMakie

    nom = NetFix.quadratic(; title = "a x^2 + b x + c = 0")
    isos = [
        isopleth(nom, s; label = "x = $(round(s.x; digits = 2))") for
        s in evaluate_all(nom; a = 1, b = -3, c = 2)
    ]
    @test length(isos) == 2
    @test length(isos[1].points) == 3          # scale, net, scale

    out = svg(nom; isopleths = isos)
    @test !occursin("NaN", out)
    @test occursin("id=\"scale-a\"", out)
    # The mesh is stroked more lightly than the scales.
    st = NomogramStyle()
    @test occursin("stroke-width=\"$(Nomograms._num(st.gridlinewidth * 4))\"", out)
    @test occursin("stroke-width=\"$(Nomograms._num(st.axislinewidth * 4))\"", out)
    @test st.gridlinewidth < st.axislinewidth

    @test nomogramplot(nom; isopleths = isos) isa Figure
end
