@testmodule ConcFix begin
    using Nomograms

    const P = 101.325
    ps(T) = 0.61094 * exp(17.625T / (T + 243.04))
    pv(W) = W * P / (622 + W)
    rh(T, W) = pv(W) / ps(T)
    enthalpy(T, W) = 1.006T + W / 1000 * (2501 + 1.86T)

    psychrometric() = concurrency_chart(
        Variable(:T; domain = (0, 45), label = "dry-bulb temperature (°C)"),
        Variable(:W; domain = (0, 30), label = "humidity ratio (g/kg)"),
        Variable(:RH; domain = (0.1, 1), ticks = 0.1:0.1:1.0) => rh,
        Variable(:h; label = "enthalpy (kJ/kg)") => enthalpy;
        region = (T, W) -> pv(W) <= ps(T) * (1 + 1e-9),
    )

    rough(Re, f) = 3.7 * (10^(-1 / (2sqrt(f))) - 2.51 / (Re * sqrt(f)))
    moody() = concurrency_chart(
        Variable(:Re, LogScale(); domain = (4e3, 1e8)),
        Variable(:f, LogScale(); domain = (0.008, 0.1)),
        Variable(:roughness; domain = (1e-6, 0.05)) => rough;
        region = (Re, f) -> rough(Re, f) >= 0,
    )

    z(u, v) = (1 + complex(u, v)) / (1 - complex(u, v))
    const RMAJ = [0.2, 0.5, 1, 2, 5, 10, 20]
    const RMIN = [0.1, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1.2, 1.4, 1.6, 1.8, 3, 4]
    smith() = concurrency_chart(
        :u => (-1, 1),
        :v => (-1, 1),
        Variable(:r; domain = (0, 50), ticks = RMAJ, minorticks = RMIN) =>
            (u, v) -> real(z(u, v)),
        Variable(
            :x;
            domain = (-50, 50),
            ticks = [0; RMAJ; -RMAJ],
            minorticks = [RMIN; -RMIN],
        ) => (u, v) -> imag(z(u, v));
        axes = false,
        region = (u, v) -> u^2 + v^2 <= 1,
        rings = [
            Variable(:angle; domain = (-180, 180), ticks = -150:30:180) =>
                θ -> (1.06cosd(θ), 1.06sind(θ)),
            Variable(:wtg; domain = (0, 0.5), ticks = 0:0.05:0.45) =>
                l -> (1.17cos(π - 4π * l), 1.17sin(π - 4π * l)),
        ],
        labelat = (r = :x => 0, x = :r => 0),
    )
end

@testitem "a reading agrees with every field's formula" tags = [:unit] setup = [ConcFix] begin
    psy = ConcFix.psychrometric()
    # From the two axes, from an axis and a field, and from two fields.
    for given in ((T = 30, W = 12), (T = 30, RH = 0.5), (RH = 0.5, h = 60), (W = 8, h = 40))
        s = evaluate(psy; given...)
        for (k, v) in pairs(given)
            @test s[k] == v
        end
        @test s.RH ≈ ConcFix.rh(s.T, s.W) rtol = 1e-9
        @test s.h ≈ ConcFix.enthalpy(s.T, s.W) rtol = 1e-9
    end
    # Air at 30 °C and half saturated carries about 13.3 g of water per kg.
    @test evaluate(psy; T = 30, RH = 0.5).W ≈ 13.28 atol = 0.01

    moody = ConcFix.moody()
    s = evaluate(moody; Re = 1e5, roughness = 1e-4)
    @test ConcFix.rough(s.Re, s.f) ≈ 1e-4 rtol = 1e-9
    @test s.f ≈ 0.0185 atol = 1e-4
end

@testitem "every crossing is a reading" tags = [:unit] begin
    c = concurrency_chart((x, y) -> (x - 0.5)^2 + y, :x => (0, 1), :y => (0, 1), :z)
    sols = evaluate_all(c; y = 0.1, z = 0.2)
    @test length(sols) == 2
    @test [s.x for s in sols] ≈ [0.5 - sqrt(0.1), 0.5 + sqrt(0.1)]
    @test_throws ArgumentError evaluate(c; y = 0.1, z = 0.2)
    # A reading can be drawn from either solution.
    r = reading(c, sols[2]; label = "right")
    @test r[:x] ≈ 0.5 + sqrt(0.1)
    @test r.point[1] ≈ c.width * (0.5 + sqrt(0.1))
end

@testitem "readings that cannot be made are refused" tags = [:unit, :validation] setup =
    [ConcFix] begin
    psy = ConcFix.psychrometric()
    msg(f) =
        try
            f()
            ""
        catch e
            e isa ArgumentError ? e.msg : rethrow()
        end
    @test occursin("two values", msg(() -> evaluate(psy; T = 30)))
    @test occursin("two values", msg(() -> evaluate(psy; T = 30, W = 10, RH = 0.5)))
    @test occursin("outside its domain", msg(() -> evaluate(psy; T = 60, RH = 0.5)))
    @test occursin("no variable :x", msg(() -> evaluate(psy; x = 1, T = 30)))
    # Beyond saturation is outside the chart's region.
    @test occursin("no point on the chart", msg(() -> evaluate(psy; T = 10, W = 25)))
    # The chart is read at a point, so an alignment line has no meaning on it.
    iso = isopleth(
        additive_nomogram(
            Variable(:a; domain = (0, 1)),
            Variable(:b; domain = (0, 1)),
            Variable(:s),
        );
        a = 0.5,
        b = 0.5,
    )
    @test occursin("readings", msg(() -> drawing(psy; isopleths = [iso])))
end

@testitem "a separable formula is drawn in straight lines" tags = [:unit] begin
    c = concurrency_chart((x, y) -> x * y, :x => (1, 10), :y => (1, 10), :z)
    # The transforms auto_nomogram recovers: logarithmic, so the axes read like log paper.
    @test c.x.variable.transform isa Nomograms.TabulatedScale
    d = drawing(c)
    zs = [sd for sd in d.scales if startswith(string(sd.name), "z_")]
    @test !isempty(zs)
    for sd in zs, p in sd.path
        a, b = sd.path[1], sd.path[end]
        chord = (b[1] - a[1], b[2] - a[2])
        off = abs(chord[1] * (p[2] - a[2]) - chord[2] * (p[1] - a[1])) / hypot(chord...)
        @test off < 1e-3
    end
    @test evaluate(c; x = 2, y = 3).z ≈ 6 rtol = 1e-9
    @test evaluate(c; x = 2, z = 12).y ≈ 6 rtol = 1e-9

    # Not separable: still a chart, with curved lines, and it still reads correctly.
    F(x, y) = sin(x) + exp(y) + x * y
    g = concurrency_chart(F, :x => (0, 2), :y => (0, 1))
    s = evaluate(g; x = 1, z = 3)
    @test F(s.x, s.y) ≈ 3 rtol = 1e-9
end

@testitem "curves thin out where they crowd" tags = [:unit] setup = [ConcFix] begin
    moody = ConcFix.moody()
    f = only(moody.fields)
    drawn(sp) = sum(
        Nomograms._length(p) for (_, _, ps) in
        Nomograms._trace_field(moody, f, NomogramStyle(; min_curve_spacing = sp)) for
        p in ps
    )
    @test drawn(0.3) > drawn(1.2) > drawn(3.0)
    # Wherever a curve is drawn, it stands at least the minimum spacing from the
    # neighbours at least as round as it is.
    st = NomogramStyle()
    gaps = Dict(v => g for (v, _, g) in Nomograms._curve_values(f.variable))
    φ = P -> Nomograms._field_at(moody, f, P)
    for (value, _, pieces) in Nomograms._trace_field(moody, f, st),
        piece in pieces,
        P in piece

        g = Nomograms._gradient(φ, P)
        @test gaps[value] / hypot(g...) >= st.min_curve_spacing * 0.98
    end
end

@testitem "labels never overlap, and the region's edge is drawn" tags = [:unit] setup =
    [ConcFix] begin
    for c in (ConcFix.psychrometric(), ConcFix.moody())
        d = drawing(
            c;
            readings = [
                reading(
                    c,
                    first(
                        evaluate_all(
                            c;
                            (
                                c.x.variable.name => sum(c.x.domain) / 2,
                                c.y.variable.name => sum(c.y.domain) / 3,
                            )...,
                        ),
                    );
                    label = "here",
                ),
            ],
        )
        st = d.style
        boxes = NTuple{4,Float64}[]
        for sd in d.scales
            isempty(sd.label) || push!(
                boxes,
                Nomograms.label_box(
                    sd.titlepos,
                    sd.label,
                    st.titlesize,
                    sd.titlealign,
                    0.0,
                ),
            )
            for t in sd.ticks
                isnothing(t.label) || push!(
                    boxes,
                    Nomograms.label_box(
                        t.labelpos,
                        t.label,
                        st.ticklabelsize,
                        t.align,
                        t.rotation,
                    ),
                )
            end
        end
        for t in d.texts
            push!(boxes, Nomograms.label_box(t.pos, t.text, st.titlesize, t.align, 0.0))
        end
        for i in eachindex(boxes), j = (i+1):lastindex(boxes)
            @test !Nomograms._overlaps(boxes[i], boxes[j], 0.0)
        end
        @test any(sd -> startswith(string(sd.name), "region_"), d.scales)
        @test length(d.isopleths) == 2
        @test occursin("here", svg(d))
    end
end

@testitem "units come back on readings" tags = [:unit] begin
    using Unitful
    bmi = concurrency_chart(
        (m, h) -> m / h^2,
        Variable(:mass; domain = (40u"kg", 150u"kg")),
        Variable(:height; domain = (1.4u"m", 2.1u"m")),
        :BMI,
    )
    s = evaluate(bmi; mass = 70u"kg", height = 175u"cm")
    @test s.BMI ≈ 70 / 1.75^2 * u"kg/m^2" rtol = 1e-9
    @test evaluate(bmi; mass = 70u"kg", BMI = 25u"kg/m^2").height ≈ sqrt(70 / 25) * u"m" rtol =
        1e-9
end

@testitem "concurrency charts render in both backends" tags = [:integration] setup =
    [ConcFix] begin
    using CairoMakie
    psy = ConcFix.psychrometric()
    r = reading(psy; T = 30, RH = 0.5, label = "30 °C, 50%")
    @test startswith(svg(psy; readings = [r]), "<?xml")
    @test nomogramplot(psy; readings = [r]) isa Figure
    @test nomogramplot(drawing(psy; readings = [r])) isa Figure
end

@testitem "a Smith chart reads impedance, angle and wavelength" tags = [:unit] setup =
    [ConcFix] begin
    smith = ConcFix.smith()
    # z = 1 + 1j has reflection coefficient (z - 1)/(z + 1) = 0.2 + 0.4j.
    s = evaluate(smith; r = 1, x = 1)
    @test s.u ≈ 0.2 atol = 1e-9
    @test s.v ≈ 0.4 atol = 1e-9
    @test s.angle ≈ atand(0.4, 0.2) atol = 1e-6
    # Wavelengths toward the generator run clockwise from 0 at the left.
    @test s.wtg ≈ (180 - atand(0.4, 0.2)) / 720 atol = 1e-6
    # A ring value and a field: found along the line out from the centre.
    back = evaluate(smith; x = 1, wtg = s.wtg)
    @test back.r ≈ 1 atol = 1e-6
    # The r = 1 circle passes through the centre, so that line meets it twice.
    @test length(evaluate_all(smith; r = 1, angle = s.angle)) == 2
    @test_throws ArgumentError evaluate(smith; angle = 30, wtg = 0.1)
end

@testitem "a chart on its own plane keeps its proportions" tags = [:unit] setup = [ConcFix] begin
    smith = ConcFix.smith()
    @test smith.height ≈ smith.width
    d = drawing(smith; readings = [reading(smith; r = 1, x = 1, label = "1 + 1j")])
    names = Set(sd.name for sd in d.scales)
    @test !(:frame_top in names) && !(:u in names) && !(:v in names)
    @test :angle in names && :wtg in names
    # Ring ticks point away from the centre all the way round.
    for n in (:angle, :wtg), t in only(sd for sd in d.scales if sd.name === n).ticks
        @test Nomograms._dot(t.normal, Nomograms._sub(t.position, smith.centre)) > 0
    end
    # A reading is read out from the centre, across the rings.
    @test d.isopleths[1].segment[1] == smith.centre
    @test length(d.isopleths) == 1
end

@testitem "labels can be set along curves, where they cross another" tags = [:unit] setup =
    [ConcFix] begin
    smith = ConcFix.smith()
    d = drawing(smith)
    inline = [t for t in d.texts if t.halo]
    @test !isempty(inline)
    @test all(t -> -pi / 2 < t.rotation <= pi / 2, inline)
    # Resistance circles are labelled on the real axis, where x = 0.
    for r in ConcFix.RMAJ
        t = findfirst(t -> t.text == Nomograms.default_tickformat(r), inline)
        isnothing(t) && continue
        @test inline[t].pos[2] ≈ smith.centre[2] atol = 0.05
    end
    @test occursin("paint-order", svg(d))
    @test_throws ArgumentError concurrency_chart(
        :u => (-1, 1),
        :v => (-1, 1),
        :w => (u, v) -> u + v;
        labelat = (q = :u => 0,),
    )
end

@testitem "a turned label is struck only by lines that cross it" tags = [:unit] begin
    pos, text, size = (0.0, 0.0), "WWWWWWWW", 2.5
    rot = pi / 4
    # The upright box of the turned text reaches this corner; the text does not.
    box = Nomograms.label_box(pos, text, size, (:center, :center), rot)
    corner = (box[3] - 0.5, box[4] - 0.5)
    @test Nomograms.segment_hits_box(corner, (box[3] + 1, box[4] + 1), box, 0.0)
    @test !Nomograms._segment_hits_label(
        corner,
        (box[3] + 1, box[4] + 1),
        pos,
        text,
        size,
        (:center, :center),
        rot,
        0.0,
    )
    @test Nomograms._segment_hits_label(
        (-1.0, 1.0),
        (1.0, -1.0),
        pos,
        text,
        size,
        (:center, :center),
        rot,
        0.0,
    )
end
