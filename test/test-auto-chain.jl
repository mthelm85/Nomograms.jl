@testitem "three inputs become a two-stage chart" tags = [:unit] begin
    vol = auto_nomogram(
        (L, W, H) -> L * W * H,
        :L => (1, 10),
        :W => (1, 10),
        :H => (1, 10),
        :V,
    )
    @test vol isa CompoundNomogram
    @test length(vol.stages) == 2
    zspan(nom, z) = (d = nom[z].domain; d[2] - d[1])
    for (L, W, H) in ((3, 4, 2), (10, 10, 10), (1.5, 7, 9), (1, 1, 1))
        @test abs(evaluate(vol; L, W, H).V - L * W * H) / zspan(vol, :V) < 1e-8
    end
    # Any direction the chain allows.
    @test evaluate(vol; L = 3, W = 4, V = 24).H ≈ 2 rtol = 1e-6
    @test evaluate(vol; L = 3, H = 2, V = 24).W ≈ 4 rtol = 1e-6

    gas = auto_nomogram(
        (n, T, V) -> n * 8.314 * T / V,
        :n => (0.5, 5),
        :T => (250, 400),
        :V => (0.01, 0.1),
        :P,
    )
    @test evaluate(gas; n = 1, T = 300, V = 0.024).P ≈ 8.314 * 300 / 0.024 rtol = 1e-6
end

@testitem "more inputs nest" tags = [:unit] begin
    four = auto_nomogram(
        (a, b, c, d) -> a * b * c * d,
        :a => (1, 5),
        :b => (1, 5),
        :c => (1, 5),
        :d => (1, 5),
        :p,
    )
    @test length(four.stages) == 3
    @test evaluate(four; a = 2, b = 3, c = 4, d = 5).p ≈ 120 rtol = 1e-6
    @test evaluate(four; a = 2, b = 3, c = 4, p = 120).d ≈ 5 rtol = 1e-6
end

@testitem "the turning lines are bare" tags = [:unit] begin
    vol = auto_nomogram(
        (L, W, H) -> L * W * H,
        :L => (1, 10),
        :W => (1, 10),
        :H => (1, 10),
        :V,
    )
    pivots = [s for s in vol.scales if startswith(string(Nomograms.name(s)), "_t")]
    @test length(pivots) == 1
    sd = only(sd for sd in drawing(vol).scales if sd.name === Nomograms.name(pivots[1]))
    @test isempty(sd.ticks)
    @test isempty(sd.label)
end

@testitem "a formula needing a net is refused, pair by pair" tags = [:unit, :validation] begin
    msg(f) =
        try
            f()
            ""
        catch e
            e isa ArgumentError ? e.msg : rethrow()
        end
    # (a + b)·c: a and b do combine, but the first stage wants their sum on a linear scale
    # and the second on a logarithmic one.
    m = msg(
        () -> auto_nomogram(
            (a, b, c) -> (a + b) * c,
            :a => (1, 10),
            :b => (1, 10),
            :c => (1, 10),
            :d,
        ),
    )
    @test occursin("cannot split d = F(a, b, c)", m)
    @test occursin(
        ":a with :b — the two stages want their combination on different scales",
        m,
    )
    @test occursin(":a with :c — the result does not depend on them through a single", m)

    # a·b + c: the same, the other way round -- logarithmic first, linear second.
    m = msg(
        () -> auto_nomogram(
            (a, b, c) -> a * b + c,
            :a => (1, 10),
            :b => (1, 10),
            :c => (1, 10),
            :z,
        ),
    )
    @test occursin("needs a net rather than a pivot", m)

    # Not monotone in one input rules out every pairing that includes it, and the reason
    # names that input.
    m = msg(
        () -> auto_nomogram(
            (a, b, c) -> (a - 5)^2 * b * c,
            :a => (0, 10),
            :b => (1, 2),
            :c => (1, 2),
            :z,
        ),
    )
    line = only(l for l in split(m, '\n') if occursin(":a with :b", l))
    @test occursin("monotone", line) && occursin(":a", line)
end

@testitem "the argument convention" tags = [:unit, :validation] begin
    # Exactly two arguments after F: two inputs, result :z.
    @test auto_nomogram((x, y) -> x + y, :x => (0, 1), :y => (0, 1)) isa Nomogram
    # Otherwise the last argument is the result -- three after F is still two inputs.
    @test variablenames(auto_nomogram((x, y) -> x + y, :x => (0, 1), :y => (0, 1), :s))[end] ==
          :s
    @test auto_nomogram(
        (x, y, w) -> x * y * w,
        :x => (1, 2),
        :y => (1, 2),
        :w => (1, 2),
        :v,
    ) isa CompoundNomogram

    @test_throws ArgumentError auto_nomogram((x,) -> x, :x => (0, 1))
    @test_throws ArgumentError auto_nomogram(
        (x, y, w) -> x * y * w,
        :x => (1, 2),
        :x => (1, 2),
        :w => (1, 2),
        :v,
    )
    @test_throws ArgumentError auto_nomogram(
        (x, y, w) -> x * y * w,
        :x => (1, 2),
        :y,
        :w => (1, 2),
        :v,
    )
end

@testitem "automatic chains render in both backends" tags = [:integration] begin
    using CairoMakie

    gas = auto_nomogram(
        (n, T, V) -> n * 8.314 * T / V,
        Variable(:n; domain = (0.5, 5), label = "amount (mol)"),
        Variable(:T; domain = (250, 400), label = "temperature (K)"),
        Variable(:V; domain = (0.01, 0.1), label = "volume (m³)"),
        Variable(:P; label = "pressure (Pa)");
        title = "Ideal gas",
    )
    isos = isopleth(gas; n = 1, T = 300, V = 0.024, label = "1 mol, 300 K, 24 L")
    @test length(isos) == 2
    out = svg(gas; isopleths = isos)
    @test !occursin("NaN", out)
    @test occursin("Ideal gas", out)
    @test nomogramplot(gas; isopleths = isos) isa Figure
end
