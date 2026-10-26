@testitem "nicestep rounds up to 1, 2, 2.5 or 5 times a power of ten" tags = [:unit] begin
    @test Nomograms.nicestep(1.0) == 1.0
    @test Nomograms.nicestep(1.4) == 2.0
    @test Nomograms.nicestep(2.3) == 2.5
    @test Nomograms.nicestep(4.0) == 5.0
    @test Nomograms.nicestep(7.0) == 10.0
    @test Nomograms.nicestep(0.013) == 0.02
    @test Nomograms.nicestep(3400.0) == 5000.0
    @test_throws ArgumentError Nomograms.nicestep(0.0)
end

@testitem "linear ticks are round numbers inside the domain" tags = [:unit] begin
    majors, minors = Nomograms.linear_ticks(0.0, 10.0)
    @test sort(majors) == collect(0.0:2.0:10.0)
    @test all(0 .<= minors .<= 10)
    @test isempty(intersect(majors, minors))
    # Subdivisions land exactly, so they print as 0.4 rather than 0.4000000000000001.
    @test 0.4 in minors

    majors, _ = Nomograms.linear_ticks(0.3, 0.7)
    @test all(0.3 .<= majors .<= 0.7)
end

@testitem "linear majors are offered roundest first" tags = [:unit] begin
    majors, _ = Nomograms.linear_ticks(0.0, 10.0)
    # 0 and 10 are multiples of the largest nice step and must be tried before 2 or 6.
    @test findfirst(==(10.0), majors) < findfirst(==(6.0), majors)
    @test findfirst(==(0.0), majors) < findfirst(==(6.0), majors)
end

@testitem "log ticks follow decades" tags = [:unit] begin
    # Every mantissa is offered however wide the range, in the order thinning should keep
    # them: decades, then 5s, then 2s, then the rest. There is no cut-off by decade count.
    majors, minors = Nomograms.log_ticks(1.0, 1e4)
    @test majors[1:5] == [1.0, 10.0, 100.0, 1000.0, 10_000.0]
    @test 20.0 in majors
    @test findfirst(==(50.0), majors) <
          findfirst(==(20.0), majors) <
          findfirst(==(30.0), majors)
    # Each mantissa is a minor too, so one whose label is dropped keeps its tick mark.
    @test 30.0 in minors && 20.0 in minors
    @test issorted(minors)

    # Over a short range too, decades first.
    majors, _ = Nomograms.log_ticks(1.0, 100.0)
    @test 1.0 in majors && 9.0 in majors && 100.0 in majors
    @test findfirst(==(10.0), majors) < findfirst(==(9.0), majors)
    @test findfirst(==(50.0), majors) < findfirst(==(30.0), majors)

    @test_throws ArgumentError Nomograms.log_ticks(0.0, 10.0)
end

@testitem "candidate ticks honour explicit values and order" tags = [:unit] begin
    v = Variable(:a; domain = (0, 10), ticks = [10.0, 0.0, 5.0, 99.0])
    majors, minors = Nomograms.candidate_ticks(v, (0.0, 10.0))
    # Out-of-domain values are dropped; the rest keep the order they were given in.
    @test majors == [10.0, 0.0, 5.0]
    @test isempty(minors)

    vlog = Variable(:a, LogScale(); domain = (1, 1000))
    majors, _ = Nomograms.candidate_ticks(vlog, (1.0, 1000.0))
    @test majors[1:4] == [1.0, 10.0, 100.0, 1000.0]
end

@testitem "log axes label what fits, from decades to every mantissa" tags = [:unit] begin
    # One candidate set, three outcomes, decided only by how much room the axis has.
    labels(nom, n) = [
        t.value for
        t in only(sd for sd in drawing(nom).scales if sd.name === n).ticks if
        !isnothing(t.label)
    ]
    mantissa(v) = round(v / 10.0^floor(log10(v) + 1e-9); digits = 6)

    # Eight decades: no room for every mantissa, but plenty for the 1-2-5 sequence. The
    # old rule counted decades and labelled the decades alone.
    wide = auto_nomogram((V, R) -> V^2 / R, :V => (1, 240), :R => (1, 1000), :P)
    ms = Set(mantissa(v) for v in labels(wide, :P))
    @test 2.0 in ms && 5.0 in ms
    @test !(3.0 in ms)

    # About two decades on a full-height axis: every mantissa fits. The old rule counted
    # floor/ceil decades -- three, for 0.127 to 12.9 -- and labelled just 1 and 10.
    two = multiplicative_nomogram(
        Variable(:s; domain = (0.5, 5)),
        Variable(:k; domain = (0.2534, 2.576)),
        Variable(:z),
    )
    zl = labels(two, :z)
    @test count(v -> 1 < v < 10, zl) >= 6
    @test count(v -> 0.1 < v < 1, zl) >= 4

    # Inside a single mantissa step there are no whole mantissas to speak of, so tenths
    # are promoted rather than leaving the axis bare between its ends.
    narrow = multiplicative_nomogram(
        Variable(:a; domain = (100, 200)),
        Variable(:b; domain = (1, 2)),
        Variable(:c),
    )
    @test 150.0 in labels(narrow, :a)
    @test length(labels(narrow, :a)) >= 6
end

@testitem "default_tickformat is compact" tags = [:unit] begin
    @test Nomograms.default_tickformat(10.0) == "10"
    @test Nomograms.default_tickformat(-3.0) == "-3"
    @test Nomograms.default_tickformat(0.5) == "0.5"
    @test Nomograms.default_tickformat(2.4) == "2.4"
    @test Nomograms.default_tickformat(1 / 3) == "0.3333"
end

@testitem "fraction labels" tags = [:unit] begin
    f = fractionformat()
    @test f.([0, 0.5, 0.375, 1.375, -1.25, 2, 3.0625]) ==
          ["0", "1/2", "3/8", "1 3/8", "-1 1/4", "2", "3 1/16"]
    # Not a whole number of sixteenths: plain decimals rather than a wrong fraction.
    @test f(0.3) == "0.3"
    @test fractionformat(16; mixed = false)(1.375) == "11/8"
    @test fractionformat(12)(1 / 3) == "1/3"
    @test_throws ArgumentError fractionformat(0)
end

@testitem "fraction candidates come coarsest first" tags = [:unit] begin
    majors, minors = Nomograms.fraction_ticks(0, 2)
    # Wholes, then halves, then quarters...
    @test majors[1:3] == [0.0, 1.0, 2.0] ||
          majors[1:3] == [2.0, 0.0, 1.0] ||
          Set(majors[1:3]) == Set([0.0, 1.0, 2.0])
    @test Set(majors[4:5]) == Set([0.5, 1.5])
    @test all(x -> x * 16 ≈ round(x * 16), majors)
    @test issorted(minors) && length(minors) == 33
    # maxden = 12 brings thirds and sixths, and no eighths.
    m12, _ = Nomograms.fraction_ticks(0, 1; maxden = 12)
    @test any(x -> x ≈ 1 / 3, m12) && any(x -> x ≈ 1 / 6, m12)
    @test !any(x -> x ≈ 1 / 8, m12)
    # A long range drops the finest fractions rather than offering thousands.
    big, _ = Nomograms.fraction_ticks(0, 1000)
    @test length(big) <= 4001
end

@testitem "integer and fraction axes label only what they allow" tags = [:unit] begin
    labels(nom, n) = [
        t.label for sd in drawing(nom).scales if sd.name === n for
        t in sd.ticks if !isnothing(t.label)
    ]
    isint(s) = occursin(r"^-?\d+$", s)
    isfrac(s) = occursin(r"^-?(\d+|\d+/\d+|\d+ \d+/\d+)$", s)

    # A count: no 2.5 teeth, and a fractional end of the domain goes unlabelled.
    gear = multiplicative_nomogram(
        Variable(:ratio; domain = (0.5, 4)),
        Variable(:driver; domain = (9.5, 40), tickstyle = :integer),
        Variable(:driven; tickstyle = :integer),
    )
    for n in (:driver, :driven)
        @test !isempty(labels(gear, n))
        @test all(isint, labels(gear, n))
    end
    @test "15" in labels(gear, :driver) && "25" in labels(gear, :driver)
    small = additive_nomogram(
        Variable(:a; domain = (0, 5), tickstyle = :integer),
        Variable(:b; domain = (0, 5)),
        Variable(:s),
    )
    @test labels(small, :a) == string.(0:5)

    # Inches: fractions, including an end like 1 5/16 that is five figures as a decimal.
    bolt = additive_nomogram(
        Variable(:grip; domain = (0.25, 1 + 5 / 16), tickstyle = :fraction),
        Variable(:nut; domain = (0.125, 0.75), tickstyle = :fraction),
        Variable(:length; tickstyle = :fraction),
    )
    for n in (:grip, :nut, :length)
        @test all(isfrac, labels(bolt, n))
    end
    @test "1 5/16" in labels(bolt, :grip)
    @test "1/2" in labels(bolt, :nut)
    # The answers are untouched: tick styles change what is drawn, not what is computed.
    @test evaluate(bolt; grip = 1.25, nut = 0.375).length ≈ 1.625
end

@testitem "tick styles survive every construction" tags = [:unit] begin
    v = Variable(:x; domain = (1, 10), tickstyle = :fraction, maxden = 8)
    @test v.maxden == 8
    w = Nomograms.rebuild(v; transform = LogScale())
    @test w.tickstyle === :fraction && w.maxden == 8 && w.tickformat(1.125) == "1 1/8"
    # auto_nomogram keeps a tick style it is given rather than choosing its own.
    a = auto_nomogram(
        (n, k) -> n * k,
        Variable(:n; domain = (1, 20), tickstyle = :integer),
        :k => (1, 3),
        :p,
    )
    @test a[:n].variable.tickstyle === :integer
    # Grid strands are labelled in the variable's own style.
    gq = grid_nomogram(
        Grid(
            Variable(:x; domain = (0.5, 2)),
            Variable(:c; domain = (0, 1), tickstyle = :fraction, maxden = 4),
            (x, c) -> (1 / (1 + x), -c / (x * (1 + x))),
        ),
        Variable(:a; domain = (0.5, 3)) => a -> (0.0, a),
        Variable(:b; domain = (-6, 0)) => b -> (1.0, b);
        relation = (x, c, a, b) -> a * x^2 + b * x + c,
    )
    strandlabels = [s.variable.label for s in gq.scales if s.kind === :grid]
    @test "1/2" in strandlabels || "1/4" in strandlabels

    @test_throws ArgumentError Variable(:x; tickstyle = :roman)
    @test_throws ArgumentError Variable(:x; tickstyle = :fraction, maxden = 0)
end

@testitem "a domain written in Rationals gets fraction ticks" tags = [:unit] begin
    style(v) = (v.tickstyle, v.maxden)
    @test style(Variable(:a; domain = (1 // 10, 10 // 10))) == (:fraction, 10)
    @test style(Variable(:b; domain = (1 // 4, 3))) == (:fraction, 4)
    @test Nomograms._fraction_denominator((1 // 3, 3 // 4)) == 12
    # All whole, but written as Rationals: still fractions, in sixteenths.
    @test style(Variable(:d; domain = (0 // 1, 3 // 1))) == (:fraction, 16)
    # A float end, plain integers, or an explicit choice: no inference.
    @test Variable(:e; domain = (1 // 3, 0.5)).tickstyle === :auto
    @test Variable(:f; domain = (1, 10)).tickstyle === :auto
    @test Variable(:g; domain = (1 // 10, 1 // 1), tickstyle = :linear).tickstyle ===
          :linear
    @test style(Variable(:h; domain = (1 // 10, 1 // 1), maxden = 20)) == (:fraction, 20)

    # Through auto_nomogram's shorthand, as a reader would write it.
    nom = auto_nomogram(
        (original, scale) -> original * scale,
        :original => (1 // 10, 10 // 10),
        :scale => (0.1, 5.0),
        :adjusted,
    )
    labels = [
        t.label for sd in drawing(nom).scales if sd.name === :original for
        t in sd.ticks if !isnothing(t.label)
    ]
    @test "3/10" in labels && "1/2" in labels
    @test !any(l -> occursin(".", l), labels)
end

@testitem "a near-whole number prints as one" tags = [:unit] begin
    f = Nomograms.default_tickformat
    @test f(10000.000000002) == "10000"
    @test f(9999.999999998) == "10000"
    @test f(-0.0) == "0"
    @test f(0.30000000000000004) == "0.3"
    @test f(2.5) == "2.5"
    @test f(1e7) == "1e7"
    # The derived end of Ohm's law's voltage axis, which used to read "1e4".
    ohm = multiplicative_nomogram(
        Variable(:I; domain = (0.001, 10)),
        Variable(:R; domain = (1, 1000)),
        Variable(:V),
    )
    labels =
        [t.label for sd in drawing(ohm).scales for t in sd.ticks if !isnothing(t.label)]
    @test !any(l -> occursin("e", l), labels)
end
