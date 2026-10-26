@testitem "a concurrent-scale nomogram adds reciprocals" tags = [:unit] begin
    nom = concurrent_nomogram(
        Variable(:R1; domain = (0, 1000)),
        Variable(:R2; domain = (0, 1000)),
        Variable(:R),
    )
    @test all(nom[:R].domain .≈ (0.0, 500.0))
    for (a, b) in ((400, 600), (100, 900), (1000, 1000), (50, 75))
        @test evaluate(nom; R1 = a, R2 = b).R ≈ a * b / (a + b) rtol = 1e-9
        # Any one of the three can be the unknown.
        r = a * b / (a + b)
        @test evaluate(nom; R1 = a, R = r).R2 ≈ b rtol = 1e-6
    end
    @test max_relation_residual(nom) < 1e-8
    # The three scales meet at their common zero.
    @test collect(point(nom[:R1], 0)) ≈ collect(point(nom[:R], 0))
    @test collect(point(nom[:R2], 0)) ≈ collect(point(nom[:R], 0))
    @test_throws ArgumentError concurrent_nomogram(
        Variable(:a; domain = (-1, 1)),
        Variable(:b; domain = (0, 1)),
        Variable(:c),
    )
end

@testitem "a proportional nomogram solves for any of its four" tags = [:unit] begin
    gas = proportional_nomogram(
        Variable(:V; domain = (1, 10)),
        Variable(:T; domain = (250, 400)),
        Variable(:nR; domain = (0.01, 0.2)),
        Variable(:p; domain = (0.5, 5)),
    )
    @test gas isa CompoundNomogram
    names = (:V, :T, :nR, :p)
    for (V, T, p) in ((5, 300, 2), (2, 260, 2), (9, 390, 4.5))
        s = evaluate(gas; V, T, p)
        @test s.nR ≈ p * V / T rtol = 1e-9
        for unknown in names
            given = NamedTuple{Tuple(n for n in names if n != unknown)}(
                Tuple(s[n] for n in names if n != unknown),
            )
            @test evaluate(gas; given...)[unknown] ≈ s[unknown] rtol = 1e-6
        end
    end
    # Two variables share each upright line, labelled on opposite sides of it.
    @test point(gas[:V], 1)[1] ≈ point(gas[:nR], 0.1)[1]
    @test point(gas[:T], 300)[1] ≈ point(gas[:p], 2)[1]
    @test gas[:V].side != gas[:nR].side && gas[:T].side != gas[:p].side
    # A reading is two lines through one point of the diagonal, drawn as one.
    isos = isopleth(gas; V = 5, T = 300, p = 2)
    @test length(isos) == 2 && isos[1].group == isos[2].group
    @test_throws ArgumentError proportional_nomogram(
        Variable(:a; domain = (1, 2)),
        Variable(:b; domain = (1, 2)),
        Variable(:c; domain = (10, 20)),
        Variable(:d; domain = (1, 2)),
    )
end
