@testitem "a variable takes its unit from its domain" tags = [:unit] begin
    using Unitful

    v = Variable(:L; domain = (1u"ft", 10u"ft"))
    @test v.unit == u"ft"
    @test v.domain == (1.0, 10.0)
    # Explicit, and mixed ends converted to the first end's unit.
    @test Variable(:L; domain = (1, 10), unit = u"ft").domain == (1.0, 10.0)
    w = Variable(:L; domain = (0.5u"m", 300u"cm"))
    @test w.unit == u"m" && w.domain == (0.5, 3.0)
    # Ticks may be quantities too.
    @test Variable(:L; domain = (1u"m", 2u"m"), ticks = [100u"cm", 1.5u"m"]).ticks ==
          [1.0, 1.5]
    @test_throws ArgumentError Variable(:L; domain = (1u"m", 2u"s"))
    # A fraction of an inch is still a fraction.
    q = Variable(:L; domain = ((1 // 2)u"inch", 3u"inch"))
    @test q.tickstyle === :fraction && q.maxden == 2 && q.unit == u"inch"
    # A unit survives every rebuild a construction makes.
    @test Nomograms.rebuild(v; transform = LogScale()).unit == u"ft"
end

@testitem "values convert in and come back in the axis's unit" tags = [:unit] begin
    using Unitful

    rect = multiplicative_nomogram(
        Variable(:w; domain = (1u"ft", 20u"ft"), label = "width"),
        Variable(:l; domain = (1u"ft", 20u"ft"), label = "length"),
        Variable(:A; unit = u"ft^2", label = "area"),
    )
    sol = evaluate(rect; w = 10u"ft", l = 12u"ft")
    @test unit(sol.A) == u"ft^2" && unit(sol.w) == u"ft"
    @test sol.A ≈ 120u"ft^2"
    # Any unit of the right dimension, and plain numbers in the axis's own unit.
    @test evaluate(rect; w = 3u"m", l = 4u"m").A ≈ uconvert(u"ft^2", 12u"m^2")
    @test evaluate(rect; w = 10, l = 12).A ≈ 120u"ft^2"
    # Inverse problems too.
    @test evaluate(rect; w = 10u"ft", A = 120u"ft^2").l ≈ 12u"ft"
    @test evaluate_all(rect; w = 10u"ft", l = 12u"ft")[1].A ≈ 120u"ft^2"

    # A solution can be handed straight back.
    iso = isopleth(rect, sol)
    @test iso[:A] ≈ 120
    @test isopleth(rect; w = 3u"m", l = 4u"m")[:w] ≈ ustrip(u"ft", 3u"m")
    @test abs(Nomograms.relation_residual(rect, sol)) < 1e-9

    # The wrong dimension, or a quantity where no unit was given, is refused by name.
    err(f) =
        try
            f()
            ""
        catch e
            sprint(showerror, e)
        end
    m = err(() -> evaluate(rect; w = 3u"s", l = 4u"ft"))
    @test occursin(":w is in ft", m)
    plain = additive_nomogram(
        Variable(:a; domain = (0, 1)),
        Variable(:b; domain = (0, 1)),
        Variable(:c),
    )
    @test occursin("`unit = u\"m\"`", err(() -> evaluate(plain; a = 1u"m", b = 0.2)))
    @test_throws ArgumentError evaluate(plain; a = "one", b = 0.2)
    # Unitless charts are untouched.
    @test evaluate(plain; a = 0.25, b = 0.5).c isa Float64
end

@testitem "axis titles carry the unit" tags = [:unit] begin
    using Unitful

    v = Variable(:L; domain = (1u"ft", 10u"ft"), label = "length")
    @test axislabel(v) == "length (ft)"
    @test axislabel(
        Variable(:L; domain = (1u"ft", 10u"ft"), label = "length", unitlabel = false),
    ) == "length"
    @test axislabel(Variable(:L; domain = (1, 10), label = "length")) == "length"
    nom = multiplicative_nomogram(
        v,
        Variable(:k; domain = (1, 2)),
        Variable(:p; unit = u"ft"),
    )
    @test "length (ft)" in [sd.label for sd in drawing(nom).scales]
    @test occursin("ft", sprint(show, MIME"text/plain"(), nom))
end

@testitem "temperatures convert through their offsets" tags = [:unit] begin
    using Unitful

    temp = additive_nomogram(
        Variable(:T; domain = (-20u"°C", 40u"°C")),
        Variable(:dT; domain = (0, 10), unit = u"K"),
        Variable(:Tout; unit = u"°C"),
    )
    sol = evaluate(temp; T = 68u"°F", dT = 5u"K")
    @test sol.T ≈ 20u"°C"
    @test sol.Tout ≈ 25u"°C"
    @test unit(sol.Tout) == u"°C"
end

@testitem "auto_nomogram checks a formula's dimensions" tags = [:unit] begin
    using Unitful

    bmi = auto_nomogram(
        (m, h) -> m / h^2,
        Variable(:mass; domain = (40u"kg", 150u"kg")),
        Variable(:height; domain = (1.4u"m", 2.1u"m")),
        Variable(:BMI; unit = u"kg/m^2"),
    )
    @test evaluate(bmi; mass = 70u"kg", height = 175u"cm").BMI ≈ 70 / 1.75^2 * u"kg/m^2" rtol =
        1e-6
    @test evaluate(bmi; mass = 154u"lb", height = 69u"inch").BMI ≈
          uconvert(u"kg/m^2", 154u"lb" / (69u"inch")^2) rtol = 1e-6

    # The result's unit comes from the formula when not given.
    @test auto_nomogram(
        (m, h) -> m / h^2,
        :m => (40u"kg", 150u"kg"),
        :h => (1.4u"m", 2.1u"m"),
        :b,
    )[:b].variable.unit == u"kg/m^2"

    # A formula with the wrong dimensions fails before any chart is built.
    m = try
        auto_nomogram(
            (m, h) -> m / h,
            Variable(:mass; domain = (40u"kg", 150u"kg")),
            Variable(:height; domain = (1.4u"m", 2.1u"m")),
            Variable(:BMI; unit = u"kg/m^2"),
        )
        ""
    catch e
        sprint(showerror, e)
    end
    @test occursin("does not come out in the units of :BMI", m)
    @test occursin("kg m^-1", m)

    # Three inputs, a constant with units of its own, a result declared in kPa.
    R = 8.314u"J/(mol*K)"
    gas = auto_nomogram(
        (n, T, V) -> n * R * T / V,
        Variable(:n; domain = (0.5u"mol", 5u"mol")),
        Variable(:T; domain = (250u"K", 400u"K")),
        Variable(:V; domain = (10u"L", 100u"L")),
        Variable(:P; unit = u"kPa"),
    )
    P = evaluate(gas; n = 1u"mol", T = 300u"K", V = 24u"L").P
    @test unit(P) == u"kPa"
    @test P ≈ uconvert(u"kPa", 1u"mol" * R * 300u"K" / 24u"L") rtol = 1e-6
    @test evaluate(gas; n = 1u"mol", T = 26.85u"°C", V = 24u"L").P ≈ P rtol = 1e-6
    @test length(isopleth(gas; n = 1u"mol", T = 300u"K", V = 24u"L")) == 2
    @test "P (kPa)" in [sd.label for sd in drawing(gas).scales]
end

@testitem "units on compound charts and nets" tags = [:unit] begin
    using Unitful

    stage1 = multiplicative_nomogram(
        Variable(:L; domain = (1u"ft", 10u"ft")),
        Variable(:W; domain = (1u"ft", 10u"ft")),
        Variable(:t; label = "", ticks = Float64[]),
    )
    stage2 = multiplicative_nomogram(
        Variable(:t; domain = stage1[:t].domain, label = "", ticks = Float64[]),
        Variable(:H; domain = (1u"ft", 10u"ft")),
        Variable(:V; unit = u"ft^3"),
    )
    box = compound(stage1, stage2)
    V = evaluate(box; L = 3u"ft", W = 4u"ft", H = 24u"inch").V
    @test V ≈ 24u"ft^3" && unit(V) == u"ft^3"
    @test length(isopleth(box; L = 1u"m", W = 4u"ft", H = 2u"ft")) == 2

    gq = net_nomogram(
        Net(
            Variable(:x; domain = (0.5u"s", 5u"s")),
            Variable(:c; domain = (0, 4)),
            (x, c) -> (1 / (1 + x), -c / (x * (1 + x))),
        ),
        Variable(:a; domain = (0.5, 3)) => a -> (0.0, a),
        Variable(:b; domain = (-6, 0)) => b -> (1.0, b);
        relation = (x, c, a, b) -> a * x^2 + b * x + c,
    )
    roots = [s.x for s in evaluate_all(gq; a = 1, b = -3, c = 2)]
    @test all(r -> unit(r) == u"s", roots)
    @test ustrip.(u"s", roots) ≈ [1, 2] rtol = 1e-6
    @test evaluate(gq; x = 1000u"ms", c = 2, a = 1).b ≈ -3 rtol = 1e-6
end
