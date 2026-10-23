@testitem "TabulatedScale interpolates and inverts" tags = [:unit] begin
    table(n) = (
        xs = collect(range(1.0, 10.0; length = n));
        Nomograms.TabulatedScale(xs, log.(xs), 1 ./ xs)
    )
    probe = range(1.0, 10.0; length = 1001)
    err(t) = maximum(abs(transform(t, x) - log(x)) for x in probe)

    # Cubic Hermite with exact derivatives is fourth order: halving the spacing cuts the
    # error about sixteenfold. Checking the order, not just one tolerance, is what says
    # the interpolation is right rather than merely fine at this size.
    @test 12 < err(table(65)) / err(table(129)) < 20
    @test err(table(1025)) < 1e-10

    t = table(65)
    for x in (1.0, 1.37, 4.2, 9.99, 10.0)
        # The inverse is the exact inverse of the interpolant, whatever its accuracy.
        @test invtransform(t, transform(t, x)) ≈ x atol = 1e-10
    end
    # Linear extrapolation just past either end, so a solver probing an ulp outside a
    # domain gets an answer rather than an exception.
    @test isfinite(transform(t, 0.999))
    @test isfinite(transform(t, 10.001))
    # Decreasing tables invert too.
    xs = collect(range(1.0, 10.0; length = 65))
    d = Nomograms.TabulatedScale(xs, -log.(xs), -1 ./ xs)
    @test invtransform(d, transform(d, 3.3)) ≈ 3.3 atol = 1e-9
end

@testitem "separability_defect tells separable from not" tags = [:unit] begin
    sd = separability_defect
    @test sd((x, y) -> x + y, (0, 10), (0, 10)) < 1e-6
    @test sd((x, y) -> x * y, (1, 100), (1, 100)) < 1e-6
    @test sd((a, b) -> 1 / (1 / a + 1 / b), (1, 10), (1, 10)) < 1e-6
    @test sd((r, n) -> (1 + r / 100)^n, (1, 20), (1, 30)) < 1e-6
    # x + y + xy = (1+x)(1+y) - 1 is separable, though it does not look it.
    @test sd((x, y) -> x + y + x * y, (0, 5), (0, 5)) < 1e-6
    # x² + xy + y is not.
    @test sd((x, y) -> x^2 + x * y + y, (0.5, 5), (0.5, 5)) > 1e-3
end

@testitem "auto_nomogram reproduces the formula" tags = [:unit] begin
    cases = [
        ((x, y) -> x + y, :x => (0, 10), :y => (0, 10), [(3, 4), (0, 10), (10, 10)]),
        ((x, y) -> x * y, :x => (1, 100), :y => (1, 100), [(3, 4), (50, 2), (100, 100)]),
        ((a, b) -> 1 / (1 / a + 1 / b), :a => (1, 10), :b => (1, 10), [(4, 6), (1, 10)]),
        ((r, n) -> (1 + r / 100)^n, :r => (1, 20), :n => (1, 30), [(5, 10), (20, 30)]),
        # Decreasing in y: the derived transforms take opposite signs and still agree.
        ((x, y) -> x / y, :x => (1, 10), :y => (1, 10), [(6, 3), (1, 10), (10, 1)]),
        ((m, h) -> m / h^2, :m => (40, 150), :h => (1.4, 2.1), [(70, 1.75), (100, 1.8)]),
    ]
    for (F, x, y, probes) in cases
        nom = auto_nomogram(F, x, y, :z)
        zspan = nom[:z].domain[2] - nom[:z].domain[1]
        for (a, b) in probes
            got = evaluate(nom; NamedTuple{(first(x), first(y))}((a, b))...).z
            @test abs(got - F(a, b)) / zspan < 1e-7
        end
        # The relation it carries is F itself, not the derived one.
        a, b = probes[1]
        sol = evaluate(nom; NamedTuple{(first(x), first(y))}((a, b))...)
        @test abs(relation_residual(nom, sol)) / zspan < 1e-7
    end
end

@testitem "auto_nomogram solves in every direction" tags = [:unit] begin
    nom = auto_nomogram((r, n) -> (1 + r / 100)^n, :r => (1, 20), :n => (1, 30), :g)
    # The rule of 72's question, read backwards off a chart built from the formula.
    @test evaluate(nom; r = 7, g = 2).n ≈ log(2) / log(1.07) rtol = 1e-7
    @test evaluate(nom; n = 10, g = 1.05^10).r ≈ 5 rtol = 1e-6
end

@testitem "it matches the chart derived by hand" tags = [:unit] begin
    auto = auto_nomogram(
        (r, n) -> (1 + r / 100)^n,
        :rate => (1, 20),
        :years => (1, 30),
        :growth,
    )
    hand = multiplicative_nomogram(
        Variable(:years; domain = (1, 30)),
        Variable(
            :rate,
            FunctionScale(r -> log(1 + r / 100), y -> 100 * expm1(y));
            domain = (1, 20),
        ),
        Variable(:growth, LogScale(ℯ)),
    )
    @test all(auto[:growth].domain .≈ hand[:growth].domain)
    for (r, n) in ((5, 10), (12, 4.5), (20, 30))
        @test evaluate(auto; rate = r, years = n).growth ≈
              evaluate(hand; rate = r, years = n).growth rtol = 1e-8
    end
end

@testitem "axes get ticks that suit their spread" tags = [:unit] begin
    style(nom, n) = Nomograms.resolved_tickstyle(nom[n].variable)
    sum = auto_nomogram((x, y) -> x + y, :x => (1, 100), :y => (1, 100))
    @test style(sum, :x) === :linear
    prod = auto_nomogram((x, y) -> x * y, :x => (1, 100), :y => (1, 100))
    @test style(prod, :x) === :log
    @test style(prod, :z) === :log
    # Over less than a decade the question does not arise.
    small = auto_nomogram((x, y) -> x * y, :x => (2, 5), :y => (2, 5))
    @test style(small, :x) === :linear
end

@testitem "what you put on a variable is kept" tags = [:unit] begin
    nom = auto_nomogram(
        (m, h) -> m / h^2,
        Variable(:mass; domain = (40, 150), label = "mass (kg)"),
        Variable(
            :height;
            domain = (1.4, 2.1),
            label = "height (m)",
            ticks = [1.5, 1.8, 2.0],
        ),
        Variable(:BMI; label = "BMI", tickstyle = :linear);
        title = "Body mass index",
    )
    @test nom[:mass].variable.label == "mass (kg)"
    @test nom[:height].variable.ticks == [1.5, 1.8, 2.0]
    @test Nomograms.resolved_tickstyle(nom[:BMI].variable) === :linear
    @test nom.title == "Body mass index"
    # ...but the transform is derived, whatever was there.
    @test nom[:mass].variable.transform isa Nomograms.TabulatedScale
end

@testitem "the result's domain is derived, or checked" tags = [:unit, :validation] begin
    nom = auto_nomogram((x, y) -> x + y, :x => (0, 10), :y => (0, 20), :s)
    @test all(nom[:s].domain .≈ (0.0, 30.0))
    @test auto_nomogram(
        (x, y) -> x + y,
        :x => (0, 10),
        :y => (0, 20),
        Variable(:s; domain = (0, 30)),
    ) isa Nomogram
    err = try
        auto_nomogram(
            (x, y) -> x + y,
            :x => (0, 10),
            :y => (0, 20),
            Variable(:s; domain = (0, 50)),
        )
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("(0, 30)", err.msg)
end

@testitem "auto_nomogram says why it cannot" tags = [:unit, :validation] begin
    msg(f) =
        try
            f()
            ""
        catch e
            e isa ArgumentError ? e.msg : rethrow()
        end
    m = msg(() -> auto_nomogram((x, y) -> x^2 + x * y + y, :x => (0.5, 5), :y => (0.5, 5)))
    @test occursin("not separable", m)
    @test occursin("determinant_nomogram", m)

    m = msg(() -> auto_nomogram((x, y) -> (x - 5)^2 + y, :x => (0, 10), :y => (0, 5)))
    @test occursin("not monotone in :x", m)

    m = msg(() -> auto_nomogram((x, y) -> x, :x => (0, 10), :y => (0, 5)))
    @test occursin("does not depend on :y", m)

    @test_throws ArgumentError auto_nomogram((x, y) -> x + y, :x, :y => (0, 1))
    # A singularity inside the domain is caught rather than drawn.
    @test_throws ArgumentError auto_nomogram((x, y) -> x / y, :x => (1, 2), :y => (-1, 1))
end

@testitem "auto charts render in both backends" tags = [:integration] begin
    using CairoMakie

    nom = auto_nomogram((a, b) -> 1 / (1 / a + 1 / b), :R1 => (1, 10), :R2 => (1, 10), :R)
    isos = [isopleth(nom; R1 = 4, R2 = 6, label = "2.4")]
    out = svg(nom; isopleths = isos)
    @test !occursin("NaN", out)
    @test nomogramplot(nom; isopleths = isos) isa Figure
end
