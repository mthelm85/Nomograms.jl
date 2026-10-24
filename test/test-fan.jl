@testmodule Fans begin
    using Nomograms

    resistors() = additive_nomogram(
        Variable(:R1, ReciprocalScale(); domain = (1, 10), label = "R₁ (Ω)"),
        Variable(:R2, ReciprocalScale(); domain = (1, 10), label = "R₂ (Ω)"),
        Variable(:R, ReciprocalScale(); label = "R (Ω)");
        title = "Resistors in parallel",
    )

    function solar()
        stage1 = multiplicative_nomogram(
            Variable(:watts; domain = (5.0, 200.0), label = "Appliance Watts"),
            Variable(:hours; domain = (0.5, 12.0), label = "Daily Runtime (h)"),
            Variable(:wh; domain = (2.5, 2400.0), ticks = Float64[], label = ""),
        )
        stage2 = multiplicative_nomogram(
            Variable(:wh; domain = (2.5, 2400.0), ticks = Float64[], label = ""),
            Variable(
                :inv_volts;
                domain = (1 / 14.0, 1 / 11.0),
                label = "System Voltage Factor",
            ),
            Variable(:ah; label = "Daily Amp-Hours (Ah)"),
        )
        return compound(stage1, stage2; title = "Solar battery budget")
    end

    # Where the line through a straight scale's two ends meets the first one's.
    function meet(a, b)
        p, q = Nomograms.point(a, a.domain[1]), Nomograms.point(a, a.domain[2])
        r, s = Nomograms.point(b, b.domain[1]), Nomograms.point(b, b.domain[2])
        d1, d2 = q .- p, s .- r
        t =
            ((r[1] - p[1]) * d2[2] - (r[2] - p[2]) * d2[1]) /
            (d1[1] * d2[2] - d1[2] * d2[1])
        return p .+ t .* d1
    end
end

@testitem "a fan gives the same answers" tags = [:unit] setup = [Fans] begin
    par = Fans.resistors()
    f = fan(par)
    # The same kind of chart comes back, just reshaped.
    @test f isa Nomogram
    @test variablenames(f) == variablenames(par)
    for a = 1:1.5:10, b = 1:1.5:10
        @test evaluate(f; R1 = a, R2 = b).R ≈ evaluate(par; R1 = a, R2 = b).R rtol = 1e-9
    end
    @test evaluate(f; R1 = 4, R = 2.4).R2 ≈ 6 rtol = 1e-9
    @test max_relation_residual(f) < 1e-9

    # A compound chart: its turning line is fanned with the rest.
    solar = Fans.solar()
    fs = fan(solar)
    @test fs isa CompoundNomogram
    for w in (5, 50, 200), t in (0.5, 4, 12), v in (1 / 14, 1 / 12, 1 / 11)
        @test evaluate(fs; watts = w, hours = t, inv_volts = v).ah ≈ w * t * v rtol = 1e-9
    end
    @test length(isopleth(fs; watts = 50, hours = 4, inv_volts = 1 / 12)) == 2

    # A net, including a reading with two answers.
    gq = net_nomogram(
        Net(
            Variable(:x; domain = (0.5, 5)),
            Variable(:c; domain = (0, 4)),
            (x, c) -> (1 / (1 + x), -c / (x * (1 + x))),
        ),
        Variable(:a; domain = (0.5, 3)) => a -> (0.0, a),
        Variable(:b; domain = (-6, 0)) => b -> (1.0, b);
        relation = (x, c, a, b) -> a * x^2 + b * x + c,
    )
    fq = fan(gq)
    @test fq isa NetNomogram
    @test [s.x for s in evaluate_all(fq; a = 1, b = -3, c = 2)] ≈ [1, 2] rtol = 1e-6
    @test max_relation_residual(fq) < 1e-6
end

@testitem "every axis becomes a spoke" tags = [:unit] setup = [Fans] begin
    for (chart, axes) in (
            (Fans.resistors(), (:R1, :R2, :R)),
            (Fans.solar(), (:watts, :wh, :ah, :hours, :inv_volts)),
        ),
        inward in (false, true)

        f = fan(chart; inward)
        centre = Fans.meet(f[axes[1]], f[axes[2]])
        for n in axes[3:end]
            @test all(isapprox.(Fans.meet(f[axes[1]], f[n]), centre; atol = 1e-6))
        end
        # The outermost spokes are `spread` apart.
        dir(n) =
            Nomograms.point(f[n], f[n].domain[2]) .- Nomograms.point(f[n], f[n].domain[1])
        angle(u, v) = acosd(clamp(sum(u .* v) / (hypot(u...) * hypot(v...)), -1, 1))
        @test maximum(angle(dir(a), dir(b)) for a in axes for b in axes) ≈ 90 atol = 1e-6
    end
end

@testitem "a fan puts the crowded end at the rim" tags = [:unit] setup = [Fans] begin
    # A reciprocal scale crowds its large values together; the map squeezes what is near
    # the hub, so those belong at the rim. The top of R1's axis is its small values.
    f = fan(Fans.resistors())
    hub(n, v) = hypot((Nomograms.point(f[n], v) .- Fans.meet(f[:R1], f[:R2]))...)
    @test hub(:R1, 1) < hub(:R1, 10)
    # Either way round can be asked for.
    out = fan(Fans.resistors(); inward = false)
    hubo(v) = hypot((Nomograms.point(out[:R1], v) .- Fans.meet(out[:R1], out[:R2]))...)
    @test hubo(1) > hubo(10)
end

@testitem "a fan's labels clear each other and every axis" tags = [:unit] setup =
    [Fans, Clash] begin
    for chart in (Fans.resistors(), Fans.solar()), spread in (60, 90, 120)
        @test isempty(Clash.clashes(drawing(fan(chart; spread))))
    end
end

@testitem "fans render in both backends" tags = [:integration] setup = [Fans] begin
    using CairoMakie

    f = fan(Fans.resistors())
    iso = isopleth(f; R1 = 4, R2 = 6, label = "2.4")
    out = svg(f; isopleths = [iso])
    @test !occursin("NaN", out)
    @test !occursin("<path", out) && !occursin("guide", out)
    @test nomogramplot(f; isopleths = [iso]) isa Figure
end

@testitem "fan argument errors" tags = [:unit, :validation] setup = [Fans] begin
    par = Fans.resistors()
    @test_throws ArgumentError fan(par; spread = 180)
    @test_throws ArgumentError fan(par; spread = 0)
    @test_throws ArgumentError fan(par; hub = 1)
    @test_throws ArgumentError fan(par; radius = -5)
    @test_throws ArgumentError fan(par; inward = :sideways)
end
