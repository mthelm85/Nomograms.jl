@testitem "align_for uses both components of the direction" tags = [:unit] begin
    af = Nomograms.align_for
    @test af((-1.0, 0.0)) == (:right, :center)
    @test af((1.0, 0.0)) == (:left, :center)
    # The case that was wrong while every axis was vertical: a label pushed straight down
    # must hang from its top edge, not straddle the point.
    @test af((0.0, -1.0)) == (:center, :top)
    @test af((0.0, 1.0)) == (:center, :bottom)
    @test af((0.7, 0.7)) == (:left, :bottom)
    @test af((-0.7, -0.7)) == (:right, :top)
end

@testitem "label_box covers the text for any alignment or rotation" tags = [:unit] begin
    lb = Nomograms.label_box
    # "ab" at size 10 is 12 wide (0.6 * size * nchars) and 10 tall.
    @test all(lb((0.0, 0.0), "ab", 10, (:left, :bottom), 0.0) .≈ (0.0, 0.0, 12.0, 10.0))
    @test all(lb((0.0, 0.0), "ab", 10, (:right, :top), 0.0) .≈ (-12.0, -10.0, 0.0, 0.0))
    @test all(lb((0.0, 0.0), "ab", 10, (:center, :center), 0.0) .≈ (-6.0, -5.0, 6.0, 5.0))
    # Quarter turn: the box turns with the text. atol, because cos(pi/2) is 6e-17 rather
    # than 0 and isapprox against an exact zero has no tolerance by default.
    @test all(
        isapprox.(
            lb((0.0, 0.0), "ab", 10, (:left, :bottom), pi / 2),
            (-10.0, 0.0, 0.0, 12.0);
            atol = 1e-12,
        ),
    )
end

@testitem "tick labels stay upright by default and turn on request" tags = [:unit] begin
    nom = additive_nomogram(
        Variable(:a; domain = (0, 10)),
        Variable(:b; domain = (0, 10)),
        Variable(:s),
    )
    flat = drawing(nom)
    @test all(t -> t.rotation == 0.0, flat.scales[1].ticks)
    @test all(t -> t.align == (:right, :center), flat.scales[1].ticks)

    turned = drawing(nom; style = NomogramStyle(; ticklabelrotation = :tangent))
    # A vertical axis read along its own direction is a quarter turn, whichever side the
    # labels are on; the two sides differ in which edge of the text meets the anchor.
    @test all(t -> t.rotation ≈ pi / 2, turned.scales[1].ticks)
    @test all(t -> t.rotation ≈ pi / 2, turned.scales[2].ticks)
    @test turned.scales[1].ticks[1].align == (:center, :bottom)
    @test turned.scales[2].ticks[1].align == (:center, :top)

    @test_throws ArgumentError drawing(
        nom;
        style = NomogramStyle(; ticklabelrotation = :sideways),
    )
end

@testitem "axis titles follow the scale off its end" tags = [:unit] begin
    # A vertical axis is unchanged: straight up from the top, centred.
    nom = additive_nomogram(
        Variable(:a; domain = (0, 10)),
        Variable(:b; domain = (0, 10)),
        Variable(:s),
    )
    sd = drawing(nom).scales[1]
    @test sd.titlepos[1] ≈ 0.0
    @test sd.titlepos[2] ≈ nom.height + NomogramStyle().titlegap
    @test sd.titlealign == (:center, :bottom)

    # A diagonal scale continues along its own direction instead of jumping to whichever
    # point happens to be highest.
    diag = Scale(Variable(:d; domain = (0, 100)), (0.0, 100.0), t -> (t, t), :right)
    flank(x, n) = Scale(Variable(n; domain = (0, 100)), (0.0, 100.0), y -> (x, y), :left)
    nom2 = Nomogram([flank(-50.0, :u), flank(150.0, :v), diag]; width = 100, height = 100)
    sd2 = drawing(nom2).scales[3]
    g = NomogramStyle().titlegap / sqrt(2)
    @test sd2.titlepos[1] ≈ 100 + g
    @test sd2.titlepos[2] ≈ 100 + g
    @test sd2.titlealign == (:left, :bottom)
end

@testmodule Crowded begin
    using Nomograms

    # Two scales close enough that their labels, which face each other across the gap,
    # cannot both be drawn.
    function chart()
        v(x, n, side) =
            Scale(Variable(n; domain = (0, 100)), (0.0, 100.0), y -> (x, y), side)
        return Nomogram(
            [v(0.0, :u, :right), v(9.0, :v, :left), v(60.0, :w, :right)];
            width = 60,
            height = 100,
        )
    end

    labelled(sd) = count(t -> !isnothing(t.label), sd.ticks)
end

@testitem "labels from different scales do not overlap" tags = [:unit] setup = [Crowded] begin
    nom = Crowded.chart()
    on = drawing(nom)
    off = drawing(nom; style = NomogramStyle(; resolve_collisions = false))

    total(d) = sum(Crowded.labelled, d.scales)
    @test total(on) < total(off)

    # Only the labels go; the tick marks themselves are untouched.
    @test [length(sd.ticks) for sd in on.scales] == [length(sd.ticks) for sd in off.scales]

    # Nothing that survived overlaps anything else that survived.
    st = on.style
    boxes = [
        Nomograms.label_box(t.labelpos, t.label, st.ticklabelsize, t.align, t.rotation)
        for sd in on.scales for t in sd.ticks if !isnothing(t.label)
    ]
    for i in eachindex(boxes), j = (i+1):lastindex(boxes)
        @test !Nomograms._overlaps(boxes[i], boxes[j], st.collisionpad)
    end
    # Turning it off really does leave overlaps, so the test above has teeth.
    offboxes = [
        Nomograms.label_box(t.labelpos, t.label, st.ticklabelsize, t.align, t.rotation)
        for sd in off.scales for t in sd.ticks if !isnothing(t.label)
    ]
    @test any(
        Nomograms._overlaps(offboxes[i], offboxes[j], st.collisionpad) for
        i in eachindex(offboxes) for j = (i+1):lastindex(offboxes)
    )
end

@testitem "axis titles outrank tick labels" tags = [:unit] setup = [Crowded] begin
    nom = Crowded.chart()
    d = drawing(nom)
    st = d.style
    titles = [
        Nomograms.label_box(sd.titlepos, sd.label, st.titlesize, sd.titlealign, 0.0) for
        sd in d.scales if !isempty(sd.label)
    ]
    @test length(titles) == 3
    for sd in d.scales, t in sd.ticks
        isnothing(t.label) && continue
        box =
            Nomograms.label_box(t.labelpos, t.label, st.ticklabelsize, t.align, t.rotation)
        @test !any(b -> Nomograms._overlaps(box, b, st.collisionpad), titles)
    end
end

@testitem "isopleth labels are nudged clear rather than dropped" tags = [:unit] begin
    nom = additive_nomogram(
        Variable(:R1, ReciprocalScale(); domain = (1, 10)),
        Variable(:R2, ReciprocalScale(); domain = (1, 10)),
        Variable(:R, ReciprocalScale());
        title = "Resistors in parallel",
    )
    # These two lines end close together: in the first release their labels overlapped.
    isos = [
        isopleth(nom; R1 = 4, R2 = 6, label = "4 || 6 = 2.4"),
        isopleth(nom; R1 = 1.5, R2 = 9, label = "1.5 || 9 = 1.29"),
    ]
    d = drawing(nom; isopleths = isos)
    st = d.style
    @test length(d.isopleths) == 2
    # Both survive: a label the caller asked for by name is never thrown away.
    @test all(id -> !isempty(id.label), d.isopleths)

    box(id) =
        Nomograms.label_box(id.labelpos, id.label, st.ticklabelsize, id.labelalign, 0.0)
    @test !Nomograms._overlaps(box(d.isopleths[1]), box(d.isopleths[2]), st.collisionpad)

    # And they clear the tick labels too.
    for sd in d.scales, t in sd.ticks
        isnothing(t.label) && continue
        tb = Nomograms.label_box(t.labelpos, t.label, st.ticklabelsize, t.align, t.rotation)
        for id in d.isopleths
            @test !Nomograms._overlaps(tb, box(id), st.collisionpad)
        end
    end
end

@testitem "the bounding box still covers rotated labels" tags = [:unit] begin
    nom = additive_nomogram(
        Variable(:a; domain = (0, 10), label = "a long axis name"),
        Variable(:b; domain = (0, 10)),
        Variable(:s);
        title = "A title",
    )
    for rot in (:horizontal, :tangent)
        d = drawing(
            nom;
            style = NomogramStyle(; ticklabelrotation = rot),
            isopleths = [isopleth(nom; a = 3, b = 4, label = "3 + 4")],
        )
        xmin, ymin, xmax, ymax = d.bbox
        st = d.style
        covers(b) =
            xmin <= b[1] + 1e-9 &&
            ymin <= b[2] + 1e-9 &&
            xmax >= b[3] - 1e-9 &&
            ymax >= b[4] - 1e-9
        for sd in d.scales
            for t in sd.ticks
                isnothing(t.label) && continue
                @test covers(
                    Nomograms.label_box(
                        t.labelpos,
                        t.label,
                        st.ticklabelsize,
                        t.align,
                        t.rotation,
                    ),
                )
            end
            @test covers(
                Nomograms.label_box(
                    sd.titlepos,
                    sd.label,
                    st.titlesize,
                    sd.titlealign,
                    0.0,
                ),
            )
        end
        for id in d.isopleths
            @test covers(
                Nomograms.label_box(
                    id.labelpos,
                    id.label,
                    st.ticklabelsize,
                    id.labelalign,
                    0.0,
                ),
            )
        end
        @test covers(
            Nomograms.label_box(d.titlepos, d.title, st.titlesize, d.titlealign, 0.0),
        )
    end
end

@testitem "both backends honour alignment and rotation" tags = [:integration] begin
    using CairoMakie

    nom = additive_nomogram(
        Variable(:a; domain = (0, 10)),
        Variable(:b; domain = (0, 10)),
        Variable(:s),
    )
    style = NomogramStyle(; ticklabelrotation = :tangent)
    out = svg(nom; style)
    # A quarter turn clockwise in SVG's y-down frame is -90 degrees.
    @test occursin("rotate(-90.00", out)
    @test occursin("text-anchor=\"middle\"", out)
    @test occursin("dominant-baseline=\"auto\"", out)
    @test occursin("dominant-baseline=\"hanging\"", out)
    @test !occursin("NaN", out)

    @test nomogramplot(nom; style) isa Figure
end
