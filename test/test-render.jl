@testmodule Render begin
    using Nomograms

    chart() = additive_nomogram(
        Variable(:R1, ReciprocalScale(); domain = (1, 10), label = "R1 (ohm)"),
        Variable(:R2, ReciprocalScale(); domain = (1, 10), label = "R2 (ohm)"),
        Variable(:R, ReciprocalScale(); label = "R (ohm)");
        title = "Resistors in parallel",
    )

    # Count opening and closing tags so a malformed document is caught without pulling
    # in an XML parser.
    function tagbalance(s, tag)
        opens = count(c -> true, eachmatch(Regex("<$tag\\b"), s))
        closes = count(c -> true, eachmatch(Regex("</$tag>"), s))
        selfclosing = count(c -> true, eachmatch(Regex("<$tag\\b[^>]*/>"), s))
        return opens - selfclosing - closes
    end
end

@testitem "svg output is a well-formed document" tags = [:unit] setup = [Render] begin
    nom = Render.chart()
    out = svg(nom; isopleths = [isopleth(nom; R1 = 4, R2 = 6, label = "4 || 6")])

    @test startswith(out, "<?xml")
    @test occursin("xmlns=\"http://www.w3.org/2000/svg\"", out)
    @test Render.tagbalance(out, "svg") == 0
    @test Render.tagbalance(out, "g") == 0
    @test Render.tagbalance(out, "text") == 0

    # The content really is there.
    @test occursin("Resistors in parallel", out)
    @test occursin("R1 (ohm)", out)
    @test occursin(">4 || 6<", out)
    @test occursin("id=\"scale-R1\"", out)
    @test occursin("id=\"isopleth-1\"", out)
    @test !occursin("NaN", out)
end

@testitem "svg escapes text and honours its options" tags = [:unit] setup = [Render] begin
    nom = additive_nomogram(
        Variable(:a; domain = (0, 10), label = "a & b < c"),
        Variable(:b; domain = (0, 10)),
        Variable(:s),
    )
    out = svg(nom)
    @test occursin("a &amp; b &lt; c", out)
    @test !occursin("a & b < c", out)

    @test occursin("fill=\"white\"", svg(nom))
    @test !occursin("<rect", svg(nom; background = "none"))

    # pxperunit scales the canvas linearly.
    w(s) = parse(Float64, match(r"width=\"([0-9.]+)\"", s)[1])
    @test w(svg(nom; pxperunit = 8)) ≈ 2 * w(svg(nom; pxperunit = 4))
end

@testitem "save_svg writes the file" tags = [:unit] setup = [Render] begin
    mktempdir() do dir
        path = joinpath(dir, "chart.svg")
        @test save_svg(path, Render.chart()) == path
        @test isfile(path)
        @test startswith(read(path, String), "<?xml")
    end
end

@testitem "nomogramplot needs Makie" tags = [:unit] begin
    # The stub only fires when the extension is absent; here CairoMakie is loaded by the
    # Makie test item, so just check the error exists for an argument nothing handles.
    @test_throws ArgumentError nomogramplot("not a nomogram")
end

@testitem "the Makie extension draws a figure" tags = [:integration] setup = [Render] begin
    using CairoMakie

    nom = Render.chart()
    iso = isopleth(nom; R1 = 4, R2 = 6, label = "4 || 6")
    fig = nomogramplot(nom; isopleths = [iso], pxperunit = 6)
    @test fig isa Figure

    ax = only(contents(fig[1, 1]))
    @test ax isa Axis
    # Three axis lines, three sets of tick segments, one isopleth line and its markers.
    @test count(p -> p isa Lines, ax.scene.plots) >= 4
    @test count(p -> p isa LineSegments, ax.scene.plots) >= 3
    @test count(p -> p isa Scatter, ax.scene.plots) >= 1

    mktempdir() do dir
        for ext in ("png", "svg", "pdf")
            path = joinpath(dir, "chart.$ext")
            save(path, fig)
            @test filesize(path) > 0
        end
    end
end

@testitem "nomogramplot! draws into an existing axis" tags = [:integration] setup = [Render] begin
    using CairoMakie

    fig = Figure()
    ax = Axis(fig[1, 1]; aspect = DataAspect())
    before = length(ax.scene.plots)
    @test nomogramplot!(ax, Render.chart()) === ax
    @test length(ax.scene.plots) > before
end
