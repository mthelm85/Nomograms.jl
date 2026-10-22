@testitem "evaluate solves in every direction" tags = [:unit] begin
    nom = additive_nomogram(
        Variable(:a; domain = (0, 10)),
        Variable(:b; domain = (0, 10)),
        Variable(:s),
    )
    # Because the solver only ever intersects a line with a curve, the unknown may be
    # any of the three.
    @test evaluate(nom; a = 3, b = 4).s ≈ 7 atol = 1e-9
    @test evaluate(nom; a = 3, s = 7).b ≈ 4 atol = 1e-9
    @test evaluate(nom; b = 4, s = 7).a ≈ 3 atol = 1e-9

    sol = evaluate(nom; a = 3, b = 4)
    @test keys(sol) == (:a, :b, :s)
    @test relation_residual(nom, sol) ≈ 0 atol = 1e-9
end

@testitem "the solver agrees with the closed form" tags = [:unit] begin
    nom = additive_nomogram(
        Variable(:R1, ReciprocalScale(); domain = (1, 10)),
        Variable(:R2, ReciprocalScale(); domain = (1, 10)),
        Variable(:R, ReciprocalScale()),
    )
    for r1 in range(1, 10; length = 11), r2 in range(1, 10; length = 11)
        @test evaluate(nom; R1 = r1, R2 = r2).R ≈ 1 / (1 / r1 + 1 / r2) rtol = 1e-9
    end
end

@testitem "evaluate argument errors" tags = [:unit, :validation] begin
    nom = additive_nomogram(
        Variable(:a; domain = (0, 10)),
        Variable(:b; domain = (0, 10)),
        Variable(:s),
    )
    @test_throws ArgumentError evaluate(nom; a = 3)                 # too few
    @test_throws ArgumentError evaluate(nom; a = 3, b = 4, s = 7)   # none unknown
    @test_throws ArgumentError evaluate(nom; a = 3, q = 4)          # no such variable
    @test_throws ArgumentError evaluate(nom; a = 3, b = 99)         # outside the domain

    # In range for its own scale, but the alignment line misses the third axis.
    nar = additive_nomogram(
        Variable(:a; domain = (0, 10)),
        Variable(:b; domain = (0, 10)),
        Variable(:s; domain = (0, 20)),
    )
    @test_throws ArgumentError evaluate(nar; a = 1, s = 19)
end

@testitem "isopleths record the crossing points" tags = [:unit] begin
    nom = additive_nomogram(
        Variable(:a; domain = (0, 10)),
        Variable(:b; domain = (0, 10)),
        Variable(:s),
    )
    iso = isopleth(nom; a = 3, b = 4, label = "3 + 4")
    @test iso.label == "3 + 4"
    @test iso[:s] ≈ 7 atol = 1e-9
    @test length(iso.points) == 3
    @test iso.points[1] == point(nom[:a], 3)

    # The three recorded points really are collinear.
    p, q, r = iso.points
    cross = (q[1] - p[1]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[1] - p[1])
    @test cross ≈ 0 atol = 1e-8
end

@testitem "align_residual changes sign across a scale" tags = [:unit] begin
    nom = additive_nomogram(
        Variable(:a; domain = (0, 10)),
        Variable(:b; domain = (0, 10)),
        Variable(:s),
    )
    p = point(nom[:a], 3)
    q = point(nom[:b], 4)
    s = nom[:s]
    @test align_residual(s, 0.0, p, q) * align_residual(s, 20.0, p, q) < 0
    @test align_residual(s, 7.0, p, q) ≈ 0 atol = 1e-8
    # Two coincident points define no line.
    @test_throws ArgumentError Nomograms.solve_on(s, p, p)
end

@testitem "a nomogram built from raw scales still solves" tags = [:unit] begin
    # No construction, no algebraic relation: three hand-placed curves, one of them bent.
    # The solver needs nothing else.
    va = Variable(:u; domain = (0, 1))
    vb = Variable(:v; domain = (0, 1))
    vc = Variable(:w; domain = (0, 1))
    scales = [
        Scale(va, (0.0, 1.0), t -> (0.0, 100t), :left),
        Scale(vb, (0.0, 1.0), t -> (100.0, 100t), :right),
        Scale(vc, (0.0, 1.0), t -> (50.0 + 10sinpi(t), 100t), :right),
    ]
    nom = Nomogram(scales; width = 100, height = 100)
    @test isnan(relation_residual(nom, (u = 0.0, v = 0.0, w = 0.0)))

    sol = evaluate(nom; u = 0.2, v = 0.8)
    p, q = point(nom[:u], 0.2), point(nom[:v], 0.8)
    @test align_residual(nom[:w], sol.w, p, q) ≈ 0 atol = 1e-8
end

@testitem "nomogram construction validates its scales" tags = [:unit, :validation] begin
    v = Variable(:a; domain = (0, 1))
    s = Scale(v, (0.0, 1.0), t -> (0.0, t), :left)
    @test_throws ArgumentError Nomogram([s, s]; width = 1, height = 1)
    @test_throws ArgumentError Nomogram([s, s, s]; width = 1, height = 1)  # duplicate names
end
