@testitem "Aqua finds nothing" tags = [:unit, :quality] begin
    using Aqua
    Aqua.test_all(Nomograms)
end

@testitem "every example runs" tags = [:integration] begin
    using CairoMakie
    dir = joinpath(pkgdir(Nomograms), "examples")
    examples = filter(endswith(".jl"), readdir(dir; join = true))
    @test !isempty(examples)
    for f in examples
        ok = try
            redirect_stdout(devnull) do
                Base.include(Module(), f)
            end
            true
        catch err
            @error "example failed" f err
            false
        end
        @test ok
    end
end

@testitem "project reshapes every kind of alignment chart" tags = [:unit] begin
    corners = [(0, 0), (70, 0), (70, 100), (0, 100)]
    target = [(14, 0), (56, 0), (70, 100), (0, 100)]

    # A chain, built by hand and found automatically. Its stages must stay aligned, or the
    # second line would not pass through the point the first one found.
    s1 = multiplicative_nomogram(
        Variable(:L; domain = (1, 10)),
        Variable(:W; domain = (1, 10)),
        Variable(:t; ticks = Float64[]),
    )
    s2 = multiplicative_nomogram(
        Variable(:t; domain = s1[:t].domain, ticks = Float64[]),
        Variable(:H; domain = (1, 10)),
        Variable(:V),
    )
    chain = compound(s1, s2)
    bent = project(chain, corners => target)
    @test bent isa CompoundNomogram
    @test evaluate(bent; L = 3, W = 4, H = 2).V ≈ 24 rtol = 1e-8
    @test evaluate(bent; L = 3, W = 4, V = 24).H ≈ 2 rtol = 1e-8
    @test collect(point(bent.stages[1][:t], 12)) ≈ collect(point(bent.stages[2][:t], 12))
    @test startswith(svg(bent; isopleths = isopleth(bent; L = 3, W = 4, H = 2)), "<?xml")

    auto = auto_nomogram(
        (a, b, c) -> a * b / c,
        :a => (1, 10),
        :b => (1, 10),
        :c => (1, 10),
        :d,
    )
    @test evaluate(project(auto, [1 0 0; 0 1 0; 0.002 0 1]); a = 2, b = 3, c = 4).d ≈ 1.5 rtol =
        1e-6

    # A grid: the surface has to move with the strands drawn from it.
    grid = grid_nomogram(
        Grid(
            Variable(:x; domain = (0.5, 5)),
            Variable(:c; domain = (0, 4)),
            (x, c) -> (1 / (1 + x), -c / (x * (1 + x))),
        ),
        Variable(:a; domain = (0.5, 3)) => a -> (0.0, a),
        Variable(:b; domain = (-6, 0)) => b -> (1.0, b);
        relation = (x, c, a, b) -> a * x^2 + b * x + c,
    )
    moved = project(grid, [1 0 0; 0 1 0; 0.002 0 1])
    @test moved isa GridNomogram
    @test [s.x for s in evaluate_all(moved; a = 1, b = -3, c = 2)] ≈ [1, 2] rtol = 1e-6
    @test max_relation_residual(moved) < 1e-6

    # The result fills the canvas it was given.
    wide = project(chain, corners => target; width = 200, height = 50)
    xs = [point(s, v)[1] for s in wide.scales for v in s.domain]
    @test maximum(xs) - minimum(xs) ≈ 200 rtol = 0.05

    @test_throws ArgumentError project(chain, zeros(3, 3))
end

@testitem "a chain's residual covers its stages" tags = [:unit] begin
    chain = auto_nomogram(
        (a, b, c) -> a * b / c,
        :a => (1, 10),
        :b => (1, 10),
        :c => (1, 10),
        :d,
    )
    @test max_relation_residual(chain) < 1e-4
end

@testitem "the package's own pivots stay out of answers" tags = [:unit] begin
    chain = auto_nomogram(
        (a, b, c) -> a * b / c,
        :a => (1, 10),
        :b => (1, 10),
        :c => (1, 10),
        :d,
    )
    s = evaluate(chain; a = 2, b = 3, c = 4)
    @test !any(n -> startswith(string(n), "_"), keys(s))
    @test Set(keys(s)) == Set((:a, :b, :c, :d))
    # And such an answer can still be drawn.
    @test length(isopleth(chain, s)) == 2

    prop = proportional_nomogram(
        Variable(:a; domain = (1, 10)),
        Variable(:b; domain = (1, 10)),
        Variable(:c; domain = (1, 10)),
        Variable(:d; domain = (1, 10)),
    )
    t = evaluate(prop; a = 2, b = 4, c = 3)
    @test Set(keys(t)) == Set((:a, :b, :c, :d))
    @test t.d ≈ 6 rtol = 1e-8
    @test length(isopleth(prop, t)) == 2

    # A pivot the caller named is theirs, and stays.
    s1 = multiplicative_nomogram(
        Variable(:L; domain = (1, 10)),
        Variable(:W; domain = (1, 10)),
        Variable(:t; ticks = Float64[]),
    )
    s2 = multiplicative_nomogram(
        Variable(:t; domain = s1[:t].domain, ticks = Float64[]),
        Variable(:H; domain = (1, 10)),
        Variable(:V),
    )
    @test evaluate(compound(s1, s2); L = 3, W = 4, H = 2).t ≈ 12 rtol = 1e-8
end

@testitem "wide and steep transforms are charted, not refused" tags = [:unit] begin
    # A fixed layout of table nodes refused all of these: each is steepest where
    # geometric nodes are widest, or spans too many decades for them.
    cases = [
        ((y, x) -> atand(y, x), (1, 100), (1, 100), (20.0, 70.0)),
        ((y, x) -> atand(y, x), (1e-3, 1e3), (1e-3, 1e3), (0.5, 40.0)),
        ((a, b) -> a * b, (1e-6, 1e6), (1e-3, 1e3), (2e3, 5e2)),
        ((a, b) -> a * b, (1e-9, 1e9), (1e-9, 1e9), (3e7, 2e5)),
        ((a, b) -> tand(a) * b, (1, 89), (1, 2), (88.0, 1.5)),
        ((a, b) -> log(a) * b, (1.001, 100), (1, 2), (1.01, 1.7)),
    ]
    for (F, xd, yd, (x, y)) in cases
        nom = auto_nomogram(F, :x => xd, :y => yd, :z)
        @test evaluate(nom; x, y).z ≈ F(x, y) rtol = 1e-4
    end
end

@testitem "mistakes are explained" tags = [:unit, :validation] begin
    msg(f) =
        try
            f()
            ""
        catch e
            e isa ArgumentError ? e.msg : "not an ArgumentError: $(typeof(e))"
        end
    A(n, d = (1, 10)) = Variable(n; domain = d)
    add = additive_nomogram(A(:a), A(:b), Variable(:c))
    disc = contour_chart((x, y) -> x + y, :x => (0, 1), :y => (0, 1), :z)

    # Using one kind of chart's tools on the other says which tool to use.
    @test occursin("isopleth", msg(() -> crossing(add; a = 2, b = 3)))
    @test occursin("crossing", msg(() -> isopleth(disc; x = 0.5, y = 0.5)))
    c = crossing(disc; x = 0.5, y = 0.5)
    @test occursin("isopleths", msg(() -> drawing(add; crossings = [c])))
    @test occursin("isopleths", msg(() -> svg(add; crossings = [c])))
    @test occursin(
        "reshapes alignment charts",
        msg(() -> project(disc, [1 0 0; 0 1 0; 0 0 1])),
    )
    @test occursin("spokes", msg(() -> fan(disc)))
    @test occursin("no separate relation", msg(() -> max_relation_residual(disc)))

    # Formulas that cannot work.
    @test occursin(
        "real number",
        msg(() -> auto_nomogram((x, y) -> complex(x, y), :x => (1, 2), :y => (1, 2))),
    )
    @test occursin(
        "real number",
        msg(() -> auto_nomogram((x, y) -> "no", :x => (1, 2), :y => (1, 2))),
    )
    @test occursin(
        "real number",
        msg(() -> contour_chart((x, y) -> [x, y], :x => (1, 2), :y => (1, 2))),
    )
    @test occursin(
        "arguments",
        msg(() -> auto_nomogram((x, y, w) -> x + y + w, :x => (1, 2), :y => (1, 2))),
    )
    @test occursin(
        "undefined",
        msg(() -> auto_nomogram((x, y) -> sqrt(x - y), :x => (0, 1), :y => (0, 1))),
    )
    # A formula's own errors are its author's to read.
    @test_throws ErrorException auto_nomogram(
        (x, y) -> error("boom"),
        :x => (1, 2),
        :y => (1, 2),
    )

    # Transforms and arguments.
    @test occursin(
        "positive",
        msg(() -> multiplicative_nomogram(A(:a, (-1, 1)), A(:b), Variable(:c))),
    )
    @test occursin(
        "range",
        msg(() -> contour_chart(:x => (0, 1), :y => (0, 1), :z => (0, 2))),
    )
    @test occursin(
        "range",
        msg(
            () ->
                contour_chart(:x => (0, 1), :y => (0, 1), :w => (x, y) -> x, :z => (0, 2)),
        ),
    )
    @test occursin(
        "true or false",
        msg(
            () -> contour_chart(
                :x => (0, 1),
                :y => (0, 1),
                :z => (x, y) -> x + y;
                region = (x, y) -> x + y,
            ),
        ),
    )
end
