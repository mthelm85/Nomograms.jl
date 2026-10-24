@testmodule Curved begin
    using Nomograms

    # Two vertical scales flanking an upward parabola: vertex at (50, 0), ends at (0, 40)
    # and (100, 40). A horizontal alignment line meets the parabola twice, which is the
    # case a sign-check at the ends of the domain gets wrong -- both ends sit on the same
    # side of the line, so the old shortcut reported "no solution" for two of them.
    function chart()
        scales = [
            Scale(Variable(:u; domain = (0, 100)), (0.0, 100.0), y -> (-50.0, y), :left),
            Scale(Variable(:v; domain = (0, 100)), (0.0, 100.0), y -> (150.0, y), :right),
            Scale(
                Variable(:w; domain = (0, 1)),
                (0.0, 1.0),
                t -> (100t, 40 * (2t - 1)^2),
                :right,
            ),
        ]
        return Nomogram(scales; width = 100, height = 100)
    end

    # 40(2t-1)^2 == 20  =>  t = (1 +/- sqrt(1/2))/2
    const TWO_CROSSINGS = ((1 - sqrt(0.5)) / 2, (1 + sqrt(0.5)) / 2)
end

@testitem "a line crossing a curved scale twice yields both solutions" tags = [:unit] setup =
    [Curved] begin
    nom = Curved.chart()
    sols = evaluate_all(nom; u = 20, v = 20)
    @test length(sols) == 2
    @test sols[1].w ≈ Curved.TWO_CROSSINGS[1] atol = 1e-8
    @test sols[2].w ≈ Curved.TWO_CROSSINGS[2] atol = 1e-8
    @test issorted([s.w for s in sols])

    # Both are genuine: each really is collinear with the two known points.
    p, q = point(nom[:u], 20), point(nom[:v], 20)
    for s in sols
        @test align_residual(nom[:w], s.w, p, q) ≈ 0 atol = 1e-6
        @test s.u == 20 && s.v == 20
    end
end

@testitem "evaluate refuses to pick one of several crossings" tags = [:unit] setup =
    [Curved] begin
    nom = Curved.chart()
    err = try
        evaluate(nom; u = 20, v = 20)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("crosses :w at 2 values", err.msg)
    @test occursin("`isopleth` with the same values", err.msg)
    @test occursin("evaluate_all", err.msg)
    # The values it cannot choose between are named, so the message is actionable.
    @test occursin("0.8536", err.msg)
end

@testitem "a single crossing still goes through evaluate" tags = [:unit] setup = [Curved] begin
    nom = Curved.chart()
    # This line leaves one of the two algebraic roots outside the scale's domain.
    sols = evaluate_all(nom; u = 100, v = 0)
    @test length(sols) == 1
    @test evaluate(nom; u = 100, v = 0).w ≈ only(sols).w
    p, q = point(nom[:u], 100), point(nom[:v], 0)
    @test align_residual(nom[:w], only(sols).w, p, q) ≈ 0 atol = 1e-6
end

@testitem "a line that misses the scale is reported as such" tags = [:unit] setup = [Curved] begin
    nom = Curved.chart()
    # y = 50 is above the parabola everywhere, so there is nothing to find.
    @test isempty(solve_all_on(nom[:w], point(nom[:u], 50), point(nom[:v], 50)))
    for f in (evaluate, evaluate_all)
        err = try
            f(nom; u = 50, v = 50)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("no solution for :w", err.msg)
    end
end

@testitem "both branches can be drawn" tags = [:unit] setup = [Curved] begin
    nom = Curved.chart()
    sols = evaluate_all(nom; u = 20, v = 20)
    isos = [isopleth(nom, s; label = "branch $i") for (i, s) in enumerate(sols)]
    @test length(isos) == 2
    @test isos[1].label == "branch 1"
    # Same alignment line, different crossing of the curved scale.
    @test isos[1].points[1] == isos[2].points[1]
    @test isos[1].points[3] != isos[2].points[3]

    d = drawing(nom; isopleths = isos)
    @test length(d.isopleths) == 2

    @test_throws ArgumentError isopleth(nom, (u = 20.0, v = 20.0))  # no :w
end

@testitem "parallel layouts are never ambiguous" tags = [:unit] begin
    # A straight scale makes the residual monotone, so there is always exactly one root
    # and evaluate_all agrees with evaluate.
    nom = additive_nomogram(
        Variable(:R1, ReciprocalScale(); domain = (1, 10)),
        Variable(:R2, ReciprocalScale(); domain = (1, 10)),
        Variable(:R, ReciprocalScale()),
    )
    for r1 in range(1, 10; length = 5), r2 in range(1, 10; length = 5)
        sols = evaluate_all(nom; R1 = r1, R2 = r2)
        @test length(sols) == 1
        @test only(sols).R ≈ evaluate(nom; R1 = r1, R2 = r2).R
    end
end

@testitem "solve_all_on still rejects a degenerate line" tags = [:unit, :validation] setup =
    [Curved] begin
    nom = Curved.chart()
    p = point(nom[:u], 20)
    @test_throws ArgumentError solve_all_on(nom[:w], p, p)
end

@testitem "a root exactly on the end of a domain is still found" tags = [:unit] begin
    # The extreme ratios of an N-chart sit precisely at the ends of the diagonal's
    # domain, and the arithmetic that derives the domain is not the arithmetic that
    # solves for the root, so the answer can land an ulp outside. Both corners must
    # still resolve.
    nom = ratio_nomogram(
        Variable(:a; domain = (1, 10)),
        Variable(:b; domain = (1, 10)),
        Variable(:c),
    )
    lo, hi = nom[:c].domain
    @test evaluate(nom; a = 1, b = 10).c ≈ lo
    @test evaluate(nom; a = 10, b = 1).c ≈ hi
    # ...and the answer is clamped onto the domain rather than left just outside it.
    @test lo <= evaluate(nom; a = 1, b = 10).c <= hi
    @test lo <= evaluate(nom; a = 10, b = 1).c <= hi

    # The same at both ends of a parallel chart.
    add = additive_nomogram(
        Variable(:x; domain = (0, 10)),
        Variable(:y; domain = (0, 10)),
        Variable(:s),
    )
    @test evaluate(add; x = 0, y = 0).s ≈ 0 atol = 1e-12
    @test evaluate(add; x = 10, y = 10).s ≈ 20
end

@testitem "isopleth draws one line through every crossing" tags = [:unit] begin
    quad = determinant_nomogram(
        Variable(:p; domain = (-6, 0)) => p -> (0.0, p),
        Variable(:q; domain = (0, 6)) => q -> (1.0, q),
        Variable(:z; domain = (0.2, 5)) => z -> (1 / (1 + z), -z^2 / (1 + z));
        relation = (p, q, z) -> z^2 + p * z + q,
    )
    # z² - 3z + 2 = 0 has both roots on the scale. evaluate still will not choose...
    @test_throws ArgumentError evaluate(quad; p = -3, q = 2)
    # ...but the line is the same for both, so isopleth draws it, marking each crossing.
    iso = isopleth(quad; p = -3, q = 2, label = "two roots")
    @test length(iso.points) == 4
    @test iso[:p] == -3 && iso[:q] == 2
    @test !haskey(iso.values, :z)
    for sol in evaluate_all(quad; p = -3, q = 2)
        @test any(
            pt -> all(isapprox.(pt, Nomograms.point(quad[:z], sol.z); atol = 1e-9)),
            iso.points,
        )
    end
    @test !occursin("NaN", svg(quad; isopleths = [iso]))
    # A single crossing is unchanged.
    one = isopleth(quad; p = -1, q = 0.25)
    @test one[:z] ≈ 0.5 rtol = 1e-6

    net = net_nomogram(
        Net(
            Variable(:x; domain = (0.5, 5)),
            Variable(:c; domain = (0, 4)),
            (x, c) -> (1 / (1 + x), -c / (x * (1 + x))),
        ),
        Variable(:a; domain = (0.5, 3)) => a -> (0.0, a),
        Variable(:b; domain = (-6, 0)) => b -> (1.0, b);
        relation = (x, c, a, b) -> a * x^2 + b * x + c,
    )
    iso = isopleth(net; a = 1, b = -3, c = 2)
    @test length(iso.points) == 4
    @test !haskey(iso.values, :x)
    @test !occursin("NaN", svg(net; isopleths = [iso]))
end
