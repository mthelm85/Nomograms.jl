@testmodule Det begin
    using Nomograms

    # The parallel construction written by hand in its own coordinates: scales at x = 0,
    # 2 and 1, the middle one at half scale, encoding u + v = w.
    parallel(; kwargs...) = determinant_nomogram(
        Variable(:u; domain = (0, 10)) => u -> (0.0, u),
        Variable(:v; domain = (0, 10)) => v -> (2.0, v),
        Variable(:w; domain = (0, 20)) => w -> (1.0, w / 2);
        relation = (u, v, w) -> u + v - w,
        kwargs...,
    )

    # z^2 + pz + q = 0. Matching (1-x3)p + x3*q = y3 against z*p + q + z^2 = 0 gives
    # x3 = 1/(1+z) and y3 = -z^2/(1+z) -- a scale no straight line can be.
    quadratic() = determinant_nomogram(
        Variable(:p; domain = (-6, 0)) => p -> (0.0, p),
        Variable(:q; domain = (0, 6)) => q -> (1.0, q),
        Variable(:z; domain = (0.2, 5)) => z -> (1 / (1 + z), -z^2 / (1 + z));
        relation = (p, q, z) -> z^2 + p * z + q,
    )
end

@testitem "curves are fitted to the canvas" tags = [:unit] setup = [Det] begin
    nom = Det.parallel()
    # The source sketch is 2 wide by 10 tall; stretching it onto 70 by 100 puts its
    # corners exactly on the canvas corners.
    @test point(nom[:u], 0) == (0.0, 0.0)
    @test point(nom[:v], 10) == (70.0, 100.0)
    @test point(nom[:w], 0)[1] ≈ 35.0

    # An affine map cannot change what a nomogram computes, so the relation survives it.
    @test evaluate(nom; u = 3, v = 4).w ≈ 7 atol = 1e-9
    @test max_relation_residual(nom) < 1e-9
end

@testitem "preserve_aspect centres instead of stretching" tags = [:unit] setup = [Det] begin
    nom = Det.parallel(; preserve_aspect = true)
    # min(70/2, 100/10) = 10, so the 2-wide sketch becomes 20 wide, centred in 70.
    @test point(nom[:u], 0)[1] ≈ 25.0
    @test point(nom[:v], 0)[1] ≈ 45.0
    @test point(nom[:u], 0)[2] ≈ 0.0
    @test point(nom[:v], 10)[2] ≈ 100.0
    # Still the same chart.
    @test evaluate(nom; u = 3, v = 4).w ≈ 7 atol = 1e-9
end

@testitem "a different canvas is the same chart" tags = [:unit] setup = [Det] begin
    a = Det.parallel()
    b = Det.parallel(; width = 300, height = 40)
    for (u, v) in ((1.0, 2.0), (7.5, 0.5), (10.0, 10.0))
        @test evaluate(a; u, v).w ≈ evaluate(b; u, v).w atol = 1e-9
    end
end

@testitem "labels face away from the middle of the chart" tags = [:unit] setup = [Det] begin
    nom = Det.parallel()
    d = drawing(nom)
    byname = Dict(sd.name => sd for sd in d.scales)
    # The outer two sit on the canvas edges, so their ticks point outwards.
    @test all(t -> t.normal[1] < 0, byname[:u].ticks)
    @test all(t -> t.normal[1] > 0, byname[:v].ticks)
    # An explicit side still wins.
    forced = determinant_nomogram(
        Variable(:u; domain = (0, 10), side = :right) => u -> (0.0, u),
        Variable(:v; domain = (0, 10)) => v -> (2.0, v),
        Variable(:w; domain = (0, 20)) => w -> (1.0, w / 2),
    )
    @test all(t -> t.normal[1] > 0, drawing(forced).scales[1].ticks)
end

@testitem "a curved scale a line meets twice" tags = [:unit] setup = [Det] begin
    nom = Det.quadratic()
    @test length(Nomograms._path(nom[:z])) > 2        # genuinely curved
    @test max_relation_residual(nom) < 1e-9

    for (p, q, roots) in ((-3, 2, (1.0, 2.0)), (-5, 4, (1.0, 4.0)), (-4, 4, (2.0,)))
        got = [s.z for s in evaluate_all(nom; p, q)]
        @test length(got) == length(roots)
        for (g, r) in zip(got, roots)
            @test g ≈ r atol = 1e-6
        end
    end
    # Two roots means evaluate declines to pick one.
    @test_throws ArgumentError evaluate(nom; p = -3, q = 2)
end

@testitem "determinant_nomogram validates its arguments" tags = [:unit, :validation] begin
    good = Variable(:a; domain = (0, 1)) => a -> (0.0, a)
    @test_throws ArgumentError determinant_nomogram(good, good)
    # Every variable needs a domain: there is no result variable to derive one from.
    @test_throws ArgumentError determinant_nomogram(
        good,
        Variable(:b) => b -> (1.0, b),
        Variable(:c; domain = (0, 1)) => c -> (0.5, c),
    )
    @test_throws ArgumentError determinant_nomogram(
        :notavariable => a -> (0.0, a),
        good,
        good,
    )
    # Scales all on one vertical line have no horizontal extent to fit.
    @test_throws ArgumentError determinant_nomogram(
        Variable(:a; domain = (0, 1)) => a -> (0.0, a),
        Variable(:b; domain = (0, 1)) => b -> (0.0, b + 1),
        Variable(:c; domain = (0, 1)) => c -> (0.0, c + 2),
    )
    # A curve that blows up inside its domain is caught rather than drawn.
    @test_throws ArgumentError determinant_nomogram(
        Variable(:a; domain = (-1, 1)) => a -> (0.0, 1 / a),
        Variable(:b; domain = (0, 1)) => b -> (1.0, b),
        Variable(:c; domain = (0, 1)) => c -> (0.5, c),
    )
end

@testitem "max_relation_residual reports rather than flatters" tags = [:unit] setup = [Det] begin
    # No relation to check against.
    plain = determinant_nomogram(
        Variable(:u; domain = (0, 10)) => u -> (0.0, u),
        Variable(:v; domain = (0, 10)) => v -> (2.0, v),
        Variable(:w; domain = (0, 20)) => w -> (1.0, w / 2),
    )
    @test isnan(max_relation_residual(plain))

    # A relation that does not match the geometry is caught.
    wrong = determinant_nomogram(
        Variable(:u; domain = (0, 10)) => u -> (0.0, u),
        Variable(:v; domain = (0, 10)) => v -> (2.0, v),
        Variable(:w; domain = (0, 20)) => w -> (1.0, w / 2);
        relation = (u, v, w) -> u * v - w,
    )
    @test max_relation_residual(wrong) > 1.0

    # And a chart where nothing is checkable says NaN rather than a clean zero.
    unreachable = determinant_nomogram(
        Variable(:u; domain = (0, 1)) => u -> (0.0, u),
        Variable(:v; domain = (0, 1)) => v -> (2.0, v),
        Variable(:w; domain = (90, 91)) => w -> (1.0, w);
        relation = (u, v, w) -> u + v - w,
    )
    @test isnan(max_relation_residual(unreachable))
end

@testitem "horizontal scales can label either side" tags = [:unit] begin
    # Ordering the two candidate normals by x alone left a horizontal scale with no way
    # to say "label above": both sides resolved to the same normal.
    up = Nomograms.side_normals((1.0, 0.0))
    @test up[1] == (0.0, -1.0)
    @test up[2] == (0.0, 1.0)
    @test up[1] != up[2]

    v = Variable(:a; domain = (0, 10))
    below = Scale(v, (0.0, 10.0), t -> (t, 50.0), :left)
    above = Scale(v, (0.0, 10.0), t -> (t, 50.0), :right)
    @test Nomograms._normal(below, 5.0)[2] < 0
    @test Nomograms._normal(above, 5.0)[2] > 0
end
