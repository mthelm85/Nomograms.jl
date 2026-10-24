@testmodule Charts begin
    using Nomograms

    simple() = additive_nomogram(
        Variable(:a; domain = (0, 10)),
        Variable(:b; domain = (0, 20)),
        Variable(:s),
    )
    resistors() = additive_nomogram(
        Variable(:R1, ReciprocalScale(); domain = (1, 10)),
        Variable(:R2, ReciprocalScale(); domain = (1, 10)),
        Variable(:R, ReciprocalScale()),
    )
    product() = multiplicative_nomogram(
        Variable(:x; domain = (1, 100)),
        Variable(:y; domain = (1, 100)),
        Variable(:z),
    )
end

@testitem "the result domain is derived from the other two" tags = [:unit] setup = [Charts] begin
    @test Charts.simple()[:s].domain == (0.0, 30.0)
    # 1/R = 1/R1 + 1/R2 over [1,10]^2 reaches R from 1/2 to 5.
    @test all(Charts.resistors()[:R].domain .≈ (0.5, 5.0))
    @test all(Charts.product()[:z].domain .≈ (1.0, 10_000.0))
end

@testitem "all three axes span the canvas height" tags = [:unit] setup = [Charts] begin
    # The middle axis is not free: its length falls out of the construction, and the
    # relation a3 * (1/a1 + 1/a2) == 1 forces it to match the outer two exactly.
    for nom in (Charts.simple(), Charts.resistors(), Charts.product())
        for s in nom.scales
            ys = [point(s, v)[2] for v in s.domain]
            @test minimum(ys) ≈ 0.0 atol = 1e-9
            @test maximum(ys) ≈ nom.height atol = 1e-9
        end
    end
end

@testitem "axis placement follows the scale factors" tags = [:unit] setup = [Charts] begin
    nom = Charts.simple()
    @test point(nom[:a], 0.0)[1] ≈ 0.0
    @test point(nom[:b], 0.0)[1] ≈ nom.width
    # a1 = H/10, a2 = H/20, so lambda = a1/(a1+a2) = 2/3.
    @test point(nom[:s], 0.0)[1] ≈ 2 / 3 * nom.width

    # Equal ranges put the result axis in the middle.
    sym = additive_nomogram(
        Variable(:a; domain = (0, 10)),
        Variable(:b; domain = (0, 10)),
        Variable(:s),
    )
    @test point(sym[:s], 0.0)[1] ≈ sym.width / 2
end

@testitem "every axis is collinear with the relation it encodes" tags = [:unit] setup =
    [Charts] begin
    # The defining property: for any valid pair of inputs, the three scale points lie on
    # one straight line. This is checked directly, without going through the solver.
    nom = Charts.resistors()
    sa, sb, sc = nom[:R1], nom[:R2], nom[:R]
    for r1 in range(1, 10; length = 7), r2 in range(1, 10; length = 7)
        r = 1 / (1 / r1 + 1 / r2)
        @test align_residual(sc, r, point(sa, r1), point(sb, r2)) ≈ 0.0 atol = 1e-8
    end
end

@testitem "an inconsistent result domain is rejected with the right one" tags = [:unit] begin
    err = try
        additive_nomogram(
            Variable(:a; domain = (0, 10)),
            Variable(:b; domain = (0, 20)),
            Variable(:s; domain = (0, 100)),
        )
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("inconsistent", err.msg)
    @test occursin("(0, 30)", err.msg)

    # The consistent one is accepted.
    @test additive_nomogram(
        Variable(:a; domain = (0, 10)),
        Variable(:b; domain = (0, 20)),
        Variable(:s; domain = (0, 30)),
    )[:s].domain == (0.0, 30.0)
end

@testitem "construction argument errors" tags = [:unit, :validation] begin
    # Only the result variable may omit its domain.
    @test_throws ArgumentError additive_nomogram(
        Variable(:a),
        Variable(:b; domain = (0, 1)),
        Variable(:s; domain = (0, 2)),
    )
    @test_throws ArgumentError additive_nomogram(
        Variable(:a; domain = (0, 1)),
        Variable(:b; domain = (0, 1)),
        Variable(:s);
        width = 0,
    )
    # A result variable with neither a domain nor an inverse cannot be placed.
    @test_throws ArgumentError additive_nomogram(
        Variable(:a; domain = (0, 1)),
        Variable(:b; domain = (0, 1)),
        Variable(:s, FunctionScale(x -> x^3 + x)),
    )
    @test_throws ArgumentError Variable(:a; domain = (5, 1))
    @test_throws ArgumentError Variable(:a; side = :middle)
end

@testitem "multiplicative charts encode a product" tags = [:unit] setup = [Charts] begin
    nom = Charts.product()
    for x in (1.5, 7.0, 42.0), y in (2.0, 9.0, 60.0)
        @test evaluate(nom; x, y).z ≈ x * y rtol = 1e-9
    end
    # Variables keep their own transforms underneath the logarithm.
    sq = multiplicative_nomogram(
        Variable(:r, PowerScale(2.0); domain = (1, 10)),
        Variable(:h; domain = (1, 10)),
        Variable(:v),
    )
    @test evaluate(sq; r = 3, h = 4).v ≈ 36.0 rtol = 1e-9
end

@testitem "the middle axis keeps clear of lopsided neighbours" tags = [:unit] begin
    # wh spans three decades, the voltage factor a tenth of one: left to fill the height
    # alike, the result axis would sit 3% of the way across.
    mk(ratio) = multiplicative_nomogram(
        Variable(:wh; domain = (2.5, 2400.0)),
        Variable(:inv_volts; domain = (1 / 14, 1 / 11)),
        Variable(:ah);
        max_modulus_ratio = ratio,
    )
    capped, full = mk(4), mk(Inf)
    x(nom, n) = Nomograms.point(nom[n], nom[n].domain[1])[1]
    W = capped.width
    @test x(full, :ah) / W < 0.05
    @test x(capped, :ah) >= W / 5 - 1e-9
    @test W - x(capped, :ah) >= W / 5 - 1e-9

    # The narrow axis is the one shortened, centred, and still on the canvas.
    y(nom, n, v) = Nomograms.point(nom[n], v)[2]
    lo, hi = y(capped, :inv_volts, 1 / 14), y(capped, :inv_volts, 1 / 11)
    @test abs(hi - lo) < capped.height / 4
    @test (lo + hi) / 2 ≈ capped.height / 2
    @test abs(y(capped, :wh, 2400.0) - y(capped, :wh, 2.5)) ≈ capped.height

    # Same chart, same answers.
    for wh in (3.0, 50.0, 2000.0), v in (1 / 13.5, 1 / 12, 1 / 11.2)
        @test evaluate(capped; wh, inv_volts = v).ah ≈ evaluate(full; wh, inv_volts = v).ah rtol =
            1e-9
        @test evaluate(capped; wh, inv_volts = v).ah ≈ wh * v rtol = 1e-9
    end
    @test max_relation_residual(capped) < 1e-9

    # Balanced ranges are untouched, and nonsense is refused.
    even = additive_nomogram(
        Variable(:a; domain = (0, 10)),
        Variable(:b; domain = (0, 20)),
        Variable(:s),
    )
    @test Nomograms.point(even[:a], 10.0)[2] ≈ even.height
    @test Nomograms.point(even[:b], 20.0)[2] ≈ even.height
    @test_throws ArgumentError additive_nomogram(
        Variable(:a; domain = (0, 10)),
        Variable(:b; domain = (0, 10)),
        Variable(:s);
        max_modulus_ratio = 0.5,
    )
end
