@testmodule Fixture begin
    using Nomograms

    chart() = additive_nomogram(
        Variable(:R1, ReciprocalScale(); domain = (1, 10), label = "R1 (ohm)"),
        Variable(:R2, ReciprocalScale(); domain = (1, 10), label = "R2 (ohm)"),
        Variable(:R, ReciprocalScale(); label = "R (ohm)");
        title = "Resistors in parallel",
    )
end

@testitem "thinning respects the minimum spacing" tags = [:unit] setup = [Fixture] begin
    nom = Fixture.chart()
    style = NomogramStyle()
    d = drawing(nom; style)
    for sd in d.scales
        majors = [t for t in sd.ticks if t.level == 1]
        for i = 2:length(majors)
            gap = hypot((majors[i].position .- majors[i-1].position)...)
            @test gap >= style.min_major_spacing - 1e-9
        end
        # Every tick sits on its own axis, at the value it claims.
        for t in sd.ticks
            @test t.position == point(nom[sd.name], t.value)
            @test hypot(t.normal...) ≈ 1
        end
    end
end

@testitem "the ends of each axis are labelled" tags = [:unit] setup = [Fixture] begin
    nom = Fixture.chart()
    d = drawing(nom)
    for sd in d.scales
        values = [t.value for t in sd.ticks]
        lo, hi = nom[sd.name].domain
        @test any(v -> isapprox(v, lo), values)
        @test any(v -> isapprox(v, hi), values)
    end
    # With an awkward domain the endpoints are not round numbers, so they are only
    # labelled because `endpoint_ticks` asks for them.
    odd = additive_nomogram(
        Variable(:a; domain = (1.3, 9.4)),
        Variable(:b; domain = (0, 10)),
        Variable(:s),
    )
    on = drawing(odd).scales[1].ticks
    off = drawing(odd; style = NomogramStyle(; endpoint_ticks = false)).scales[1].ticks
    @test any(t -> isapprox(t.value, 1.3), on)
    @test any(t -> isapprox(t.value, 9.4), on)
    @test !any(t -> isapprox(t.value, 1.3), off)
    @test !any(t -> isapprox(t.value, 9.4), off)
end

@testitem "labels go on the requested side" tags = [:unit] setup = [Fixture] begin
    nom = Fixture.chart()
    d = drawing(nom)
    byname = Dict(sd.name => sd for sd in d.scales)
    # Outer axes label outwards; the parallel-axis construction sets that up.
    @test all(t -> t.normal[1] < 0, byname[:R1].ticks)
    @test all(t -> t.normal[1] > 0, byname[:R2].ticks)

    flipped = additive_nomogram(
        Variable(:a; domain = (0, 10), side = :right),
        Variable(:b; domain = (0, 10), side = :left),
        Variable(:s),
    )
    df = drawing(flipped)
    @test all(t -> t.normal[1] > 0, df.scales[1].ticks)
    @test all(t -> t.normal[1] < 0, df.scales[2].ticks)
end

@testitem "straight axes collapse to a chord, curved ones do not" tags = [:unit] setup =
    [Fixture] begin
    d = drawing(Fixture.chart())
    @test all(sd -> length(sd.path) == 2, d.scales)

    curved = Scale(
        Variable(:w; domain = (0, 1)),
        (0.0, 1.0),
        t -> (50 + 10sinpi(t), 100t),
        :right,
    )
    @test length(Nomograms._path(curved)) > 2
end

@testitem "the bounding box contains everything drawn" tags = [:unit] setup = [Fixture] begin
    nom = Fixture.chart()
    iso = isopleth(nom; R1 = 4, R2 = 6, label = "4 || 6")
    d = drawing(nom; isopleths = [iso])
    xmin, ymin, xmax, ymax = d.bbox
    @test xmin < xmax && ymin < ymax

    inside(p) = xmin - 1e-9 <= p[1] <= xmax + 1e-9 && ymin - 1e-9 <= p[2] <= ymax + 1e-9
    for sd in d.scales
        @test all(inside, sd.path)
        @test all(t -> inside(t.position), sd.ticks)
        @test inside(sd.titlepos)
    end
    for id in d.isopleths
        @test inside(id.segment[1]) && inside(id.segment[2])
        @test all(inside, id.points)
    end
    @test inside(d.titlepos)
    # Axis titles sit above the axes, and the chart title above those.
    @test d.titlepos[2] > maximum(sd.titlepos[2] for sd in d.scales)
end

@testitem "isopleth geometry overhangs the outer scales" tags = [:unit] setup = [Fixture] begin
    nom = Fixture.chart()
    iso = isopleth(nom; R1 = 4, R2 = 6)
    style = NomogramStyle()
    d = drawing(nom; style, isopleths = [iso])
    id = only(d.isopleths)
    @test id.segment[1][1] < minimum(p[1] for p in id.points)
    @test id.segment[2][1] > maximum(p[1] for p in id.points)
    @test length(id.points) == 3
end

@testitem "a derived end is left unlabelled unless it is round" tags = [:unit] begin
    nom = multiplicative_nomogram(
        Variable(:a; domain = (1.5, 10)),
        Variable(:b; domain = (1.11, 9)),
        Variable(:c),
    )
    @test nom[:c].domain[1] ≈ 1.665
    ticks = Dict(sd.name => [t.value for t in sd.ticks] for sd in drawing(nom).scales)
    # Declared ends keep their ticks, up to three significant figures.
    @test all(v -> any(t -> t ≈ v, ticks[:a]), (1.5, 10))
    @test all(v -> any(t -> t ≈ v, ticks[:b]), (1.11, 9))
    # The derived 1.665 is wherever the arithmetic landed; 90 is a round number.
    @test !any(t -> t ≈ 1.665, ticks[:c])
    @test any(t -> t ≈ 90, ticks[:c])

    @test Nomograms._readable_end(0.0) && Nomograms._readable_end(-250)
    @test Nomograms._readable_end(1.01) && !Nomograms._readable_end(237.4)
    @test !Nomograms._readable_end(1.883e6)
end

@testitem "ticks are graded by importance" tags = [:unit] begin
    len(d, n) =
        Dict(t.value => t.length for sd in d.scales if sd.name === n for t in sd.ticks)
    st = NomogramStyle()
    bolt = additive_nomogram(
        Variable(:grip; domain = (0.25, 3), tickstyle = :fraction),
        Variable(:nut; domain = (0.125, 0.75), tickstyle = :fraction),
        Variable(:length; tickstyle = :fraction),
    )
    g = len(drawing(bolt), :grip)
    # Wholes, halves, quarters, eighths, sixteenths: each shorter than the last.
    @test g[2.0] > g[1.5] > g[1.25] > g[1.125] > g[1.0625]
    @test g[2.0] ≈ st.majorticklength && g[1.0625] ≈ st.minorticklength
    # Axis ends are drawn in full whatever their tier.
    @test g[0.25] ≈ st.majorticklength
    # Labels stay in one column however their ticks are graded.
    d = drawing(bolt)
    grip = only(sd for sd in d.scales if sd.name === :grip)
    xs = [t.labelpos[1] for t in grip.ticks if !isnothing(t.label)]
    @test maximum(xs) - minimum(xs) < 1e-9

    logs = multiplicative_nomogram(
        Variable(:a; domain = (1, 100)),
        Variable(:b; domain = (1, 10)),
        Variable(:c),
    )
    l = len(drawing(logs), :a)
    @test l[10.0] > l[20.0] > minimum(values(l))
    @test minimum(values(l)) ≈ st.minorticklength

    # graded_ticks = false gives the two fixed lengths back.
    flat = len(drawing(bolt; style = NomogramStyle(; graded_ticks = false)), :grip)
    @test Set(round.(values(flat); digits = 9)) ⊆
          Set(round.([st.majorticklength, st.minorticklength]; digits = 9))
end
