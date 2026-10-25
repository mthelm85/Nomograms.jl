@testitem "auto_nomogram puts the roots of a quadratic on a curve" tags = [:unit] begin
    root(p, q) = (-p + sqrt(p^2 - 4q)) / 2
    nom = auto_nomogram(root, :p => (-7, -5), :q => (1, 6), :root)
    # The result is the curved scale; the inputs stay straight.
    @test length(Nomograms._path(nom[:root])) > 2
    @test length(Nomograms._path(nom[:p])) == 2
    @test length(Nomograms._path(nom[:q])) == 2
    for p in (-6.8, -6.1, -5.3), q in (1.5, 3.2, 5.7)
        got = Nomograms._evaluate(nom; p, q)[:root]
        @test isapprox(got, root(p, q); rtol = 1e-5)
    end
end

@testitem "auto_nomogram finds which variable goes on the curve" tags = [:unit] begin
    # Linear in y and z with coefficients in x: x is the curved one.
    F(x, y) = x^2 + x * y + y
    nom = auto_nomogram(F, :x => (0.5, 3), :y => (0.5, 3), :z)
    @test length(Nomograms._path(nom[:x])) > 2
    @test length(Nomograms._path(nom[:y])) == 2
    @test length(Nomograms._path(nom[:z])) == 2
    for x in (0.6, 1.7, 2.9), y in (0.6, 1.7, 2.9)
        @test isapprox(Nomograms._evaluate(nom; x, y)[:z], F(x, y); rtol = 1e-5)
    end
end

@testitem "separable formulas still get three straight scales" tags = [:unit] begin
    nom = auto_nomogram((x, y) -> x * y, :x => (1, 10), :y => (1, 10), :z)
    @test all(n -> length(Nomograms._path(nom[n])) == 2, (:x, :y, :z))
end

@testitem "a formula that cannot be differentiated is explained" tags = [:unit] begin
    f(x::Float64, y::Float64) = x^2 + x * y + y
    m = try
        auto_nomogram(f, :x => (0.5, 3), :y => (0.5, 3))
        ""
    catch e
        e isa ArgumentError ? e.msg : rethrow()
    end
    @test occursin("could not be differentiated automatically", m)
end

@testitem "a curved chart's declared result domain is checked" tags = [:unit, :validation] begin
    root(p, q) = (-p + sqrt(p^2 - 4q)) / 2
    err = try
        auto_nomogram(root, :p => (-7, -5), :q => (1, 6), Variable(:root; domain = (0, 10)))
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("does not match", err.msg)
    nom = auto_nomogram(
        root,
        :p => (-7, -5),
        :q => (1, 6),
        Variable(:root; domain = (3, (7 + sqrt(45)) / 2)),
    )
    @test nom isa Nomogram
end
