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
    # "ab" at size 10 is as wide as its glyphs and 10 tall.
    w = Nomograms._textwidth("ab", 10)
    @test w ≈ 2 * 5.562 atol = 0.01
    @test all(lb((0.0, 0.0), "ab", 10, (:left, :bottom), 0.0) .≈ (0.0, 0.0, w, 10.0))
    @test all(lb((0.0, 0.0), "ab", 10, (:right, :top), 0.0) .≈ (-w, -10.0, 0.0, 0.0))
    @test all(
        lb((0.0, 0.0), "ab", 10, (:center, :center), 0.0) .≈ (-w / 2, -5.0, w / 2, 5.0),
    )
    # Quarter turn: the box turns with the text. atol, because cos(pi/2) is 6e-17 rather
    # than 0 and isapprox against an exact zero has no tolerance by default.
    @test all(
        isapprox.(
            lb((0.0, 0.0), "ab", 10, (:left, :bottom), pi / 2),
            (-10.0, 0.0, 0.0, w);
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

@testitem "segment_hits_box clips correctly" tags = [:unit] begin
    hits = Nomograms.segment_hits_box
    box = (0.0, 0.0, 10.0, 10.0)
    @test hits((-5.0, 5.0), (15.0, 5.0), box, 0)        # straight through
    @test hits((5.0, 5.0), (20.0, 5.0), box, 0)         # starting inside
    @test hits((-5.0, -5.0), (15.0, 15.0), box, 0)      # diagonal
    @test !hits((-5.0, 20.0), (15.0, 20.0), box, 0)     # parallel, clear above
    @test !hits((-5.0, 5.0), (-1.0, 5.0), box, 0)       # stops short
    @test !hits((20.0, 0.0), (20.0, 10.0), box, 0)      # parallel, clear to the right
    # Corner cases: just missing, and brought into contact by the pad.
    @test !hits((11.0, -5.0), (11.0, 15.0), box, 0)
    @test hits((11.0, -5.0), (11.0, 15.0), box, 2)
end

@testitem "a label an alignment line runs through is dropped" tags = [:unit] begin
    quad = determinant_nomogram(
        Variable(:p; domain = (-6, 0)) => p -> (0.0, p),
        Variable(:q; domain = (0, 6)) => q -> (1.0, q),
        Variable(:z; domain = (0.2, 5)) => z -> (1 / (1 + z), -z^2 / (1 + z)),
    )
    isos = [isopleth(quad, s) for s in evaluate_all(quad; p = -3, q = 2)]
    zlabels(st) = [
        t.value for t in only(
            sd for
            sd in drawing(quad; style = st, isopleths = isos).scales if sd.name === :z
        ).ticks if !isnothing(t.label)
    ]

    on = zlabels(NomogramStyle())
    off = zlabels(NomogramStyle(; avoid_isopleths = false))
    @test 3.0 in off
    @test !(3.0 in on)
    @test length(on) == length(off) - 1

    # The tick mark stays; only its number goes.
    d = drawing(quad; isopleths = isos)
    zd = only(sd for sd in d.scales if sd.name === :z)
    @test any(t -> t.value == 3.0, zd.ticks)

    # With no line to run through anything, nothing is lost.
    bare = [
        t.value for
        t in only(sd for sd in drawing(quad).scales if sd.name === :z).ticks if
        !isnothing(t.label)
    ]
    @test 3.0 in bare
    @test bare == off
end

@testitem "the overhang does not cost labels" tags = [:unit] begin
    # The line is drawn past the outermost scales as a flourish, and that flourish falls
    # exactly where those scales' labels are -- including the label for the value being
    # read. Only the span between the crossings counts.
    nom = additive_nomogram(
        Variable(:a; domain = (0, 10)),
        Variable(:b; domain = (0, 10)),
        Variable(:s),
    )
    d = drawing(nom; isopleths = [isopleth(nom; a = 4, b = 6)])
    kept(n) = [
        t.value for
        t in only(sd for sd in d.scales if sd.name === n).ticks if !isnothing(t.label)
    ]
    @test 4.0 in kept(:a)
    @test 6.0 in kept(:b)
end

@testmodule Clash begin
    using Nomograms

    # Nothing in a drawing may overlap: titles against titles, any text against another
    # scale's axis line, tick labels against titles. Returns the offending pairs.
    function clashes(d)
        st = d.style
        pad = st.collisionpad
        box(pos, text, size, align, rot = 0.0) =
            Nomograms.label_box(pos, text, size, align, rot)
        texts = Tuple{Int,NTuple{4,Float64},String}[]
        for (j, sd) in enumerate(d.scales)
            isempty(sd.label) || push!(
                texts,
                (j, box(sd.titlepos, sd.label, st.titlesize, sd.titlealign), sd.label),
            )
            for t in sd.ticks
                isnothing(t.label) && continue
                push!(
                    texts,
                    (
                        j,
                        box(t.labelpos, t.label, st.ticklabelsize, t.align, t.rotation),
                        t.label,
                    ),
                )
            end
        end
        bad = String[]
        for a = 1:length(texts), b = (a+1):length(texts)
            Nomograms._overlaps(texts[a][2], texts[b][2], pad) &&
                push!(bad, "$(texts[a][3]) / $(texts[b][3])")
        end
        for (j, bx, txt) in texts, (k, sd) in enumerate(d.scales)
            (k == j || sd.kind !== :scale) && continue
            for i = 1:(length(sd.path)-1)
                if Nomograms.segment_hits_box(sd.path[i], sd.path[i+1], bx, pad)
                    push!(bad, "$txt across axis $(sd.name)")
                    break
                end
            end
        end
        return bad
    end
end

@testitem "titles never overlap, however close the axes" tags = [:unit] setup = [Clash] begin
    # The chart this was reported on: a compound whose second stage puts its result 2.4
    # units from the pivot line, with titles a good deal wider than that.
    stage1 = multiplicative_nomogram(
        Variable(:watts; domain = (5.0, 200.0), label = "Appliance Watts"),
        Variable(:hours; domain = (0.5, 12.0), label = "Daily Runtime (h)"),
        Variable(
            :wh;
            domain = (2.5, 2400.0),
            ticks = Float64[],
            label = "Turning Line (Watt-Hours)",
        ),
    )
    stage2 = multiplicative_nomogram(
        Variable(
            :wh;
            domain = (2.5, 2400.0),
            ticks = Float64[],
            label = "Turning Line (Watt-Hours)",
        ),
        Variable(
            :inv_volts;
            domain = (1 / 14.0, 1 / 11.0),
            label = "System Voltage Factor",
        ),
        Variable(:ah; label = "Daily Amp-Hours (Ah)"),
    )
    solar = compound(stage1, stage2)
    isos =
        isopleth(solar; watts = 50.0, hours = 4.0, inv_volts = 1 / 12.0, label = "16.7 Ah")
    d = drawing(solar; isopleths = isos)
    @test isempty(Clash.clashes(d))
    # Every title still sits at an end of its own axis rather than drifting sideways.
    for sd in d.scales
        isempty(sd.label) && continue
        @test abs(sd.titlepos[1] - sd.path[1][1]) < 1e-9
    end

    # Worse still: five long titles on axes 3 units apart, which no single row can hold.
    v(x, n) = Scale(
        Variable(n; domain = (0, 100), label = "A rather long title for $(n)"),
        (0.0, 100.0),
        y -> (x, y),
        :right,
    )
    packed = Nomogram([v(3.0 * k, Symbol(:s, k)) for k = 0:4]; width = 12, height = 100)
    @test isempty(Clash.clashes(drawing(packed)))

    # Curved and diagonal scales, whose titles leave along their own direction.
    quad = determinant_nomogram(
        Variable(:p; domain = (-6, 0), label = "the linear coefficient") =>
            p -> (0.0, p),
        Variable(:q; domain = (0, 6), label = "the constant term") => q -> (1.0, q),
        Variable(:z; domain = (0.2, 5), label = "a root of the quadratic") =>
            z -> (1 / (1 + z), -z^2 / (1 + z));
        relation = (p, q, z) -> z^2 + p * z + q,
    )
    @test isempty(Clash.clashes(drawing(quad)))
    speed = ratio_nomogram(
        Variable(:distance; domain = (1, 10), label = "distance travelled (km)"),
        Variable(:time; domain = (1, 10), label = "time taken (h)"),
        Variable(:speed; label = "average speed (km/h)"),
    )
    @test isempty(Clash.clashes(drawing(speed)))
end

@testitem "nudged labels never settle on an overlap" tags = [:unit] begin
    # Far more labels wanting one spot than a few nudges can separate.
    st = NomogramStyle()
    occupied = NTuple{4,Float64}[]
    for k = 1:40
        Nomograms.nudge_clear(
            "label $k",
            (0.0, 0.0),
            (0.0, 1.0),
            st.ticklabelsize,
            (:left, :center),
            st,
            occupied,
        )
    end
    @test length(occupied) == 40
    @test !any(
        Nomograms._overlaps(occupied[i], occupied[j], st.collisionpad) for
        i in eachindex(occupied) for j = (i+1):lastindex(occupied)
    )
end

@testitem "label widths follow the font" tags = [:unit] begin
    w(s) = Nomograms._textwidth(s, 10.0)
    # Proportional, as Helvetica is: capitals and wide letters take more room than digits,
    # and narrow letters far less.
    @test w("WWWW") > 1.5 * w("0000")
    @test w("iiii") < 0.5 * w("0000")
    @test w("WAVELENGTH") > w("wavelength")
    # Digits are 0.556 em; a fraction label is narrower than a flat estimate.
    @test w("1000") ≈ 4 * 5.562 atol = 0.01
    @test w("1 3/8") < 0.6 * 10 * 5
    # Characters outside the table fall back, wide ones to a full em.
    @test w("ж") ≈ 6.0
    @test w("漢") ≈ 10.0
end

@testitem "capital titles do not overlap" tags = [:unit] setup = [Clash] begin
    # Titles heavy in wide capitals, on close axes: the estimate this replaced put these
    # boxes too narrow and let them collide.
    nom = additive_nomogram(
        Variable(:a; domain = (0, 10), label = "MAXIMUM WAVELENGTH"),
        Variable(:b; domain = (0, 10), label = "WAVE HEIGHT (M)"),
        Variable(:c; label = "WINDOW WIDTH (MM)");
        width = 40,
    )
    @test isempty(Clash.clashes(drawing(nom)))
end
