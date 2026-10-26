"""
    DEFAULT_ISOPLETH_COLORS

Blue, orange and aqua, in that order. All the isopleths on a chart are visible at once, so
every pair of colours has to stay distinguishable under colour-vision deficiency, and no
fourth hue does once every pair is taken into account. The colours were checked with a
validator: the worst pair is ΔE 9.2 under deuteranopia and 24.0 in normal vision, against a
target of at least 8.

d3's `category10` palette is not used because its orange and green are only ΔE 0.7 apart
under protanopia, which makes them look the same.
"""
const DEFAULT_ISOPLETH_COLORS = ("#2a78d6", "#eb6834", "#1baf7a")

"""
    DEFAULT_ISOPLETH_DASHES

Dash patterns, cycled one step behind the colours, so a fourth isopleth is dashed blue
instead of a second plain blue. Lines are never told apart by colour alone, so a chart can
carry more lines than there are safe colours.
"""
const DEFAULT_ISOPLETH_DASHES = (:solid, :dash, :dot, :dashdot)

"""
    isopleth_style(i, colors, dashes) -> (color, dash)

Colour and dash for the `i`-th isopleth: the colour advances every line, the dash every
full cycle of colours.
"""
function isopleth_style(i::Integer, colors, dashes)
    nc = length(colors)
    return (colors[mod1(i, nc)], dashes[mod1(cld(i, nc), length(dashes))])
end

# Dash lengths are given in chart units so they scale with the drawing, like everything
# else in NomogramStyle.
function _dasharray(dash::Symbol, k::Real)
    dash === :solid && return ""
    segs = if dash === :dash
        (1.8, 1.1)
    elseif dash === :dot
        (0.25, 0.9)
    elseif dash === :dashdot
        (1.8, 0.9, 0.25, 0.9)
    else
        throw(ArgumentError("unknown dash pattern :$dash"))
    end
    return """ stroke-dasharray="$(join((_num(s * k) for s in segs), " "))\""""
end

function _esc(s::AbstractString)
    out = replace(s, '&' => "&amp;", '<' => "&lt;", '>' => "&gt;")
    return replace(out, '"' => "&quot;")
end

_num(x::Real) = @sprintf("%.2f", x)

_svg_anchor(h::Symbol) = h === :left ? "start" : h === :right ? "end" : "middle"
_svg_baseline(v::Symbol) = v === :bottom ? "auto" : v === :top ? "hanging" : "central"

# A text element placed from the layout's anchor point, alignment and rotation. The
# rotation is negated because SVG's y axis runs downwards while the layout's runs up.
function _text(
    io::IO,
    indent::AbstractString,
    X,
    Y,
    pos::Point,
    text::AbstractString,
    size::Real,
    align::NTuple{2,Symbol},
    rotation::Real,
    fill::AbstractString,
    weight::AbstractString = "",
    extra::AbstractString = "",
)
    x, y = X(pos[1]), Y(pos[2])
    transform = if iszero(rotation)
        ""
    else
        """ transform="rotate($(_num(-rad2deg(rotation))) $(_num(x)) $(_num(y)))\""""
    end
    bold = isempty(weight) ? "" : """ font-weight="$weight\""""
    println(
        io,
        """$indent<text x="$(_num(x))" y="$(_num(y))" font-size="$(_num(size))" """ *
        """text-anchor="$(_svg_anchor(align[1]))" """ *
        """dominant-baseline="$(_svg_baseline(align[2]))" """ *
        """fill="$(_esc(fill))"$bold$transform$extra>$(_esc(text))</text>""",
    )
    return nothing
end

"""
    svg(nom::AbstractNomogram; kwargs...) -> String
    svg(d::NomogramDrawing; kwargs...) -> String

Render a nomogram as a standalone SVG document, returned as a `String`.

This backend needs nothing beyond the standard library, so vector output works without any
other packages. To put the chart in a Makie `Figure`, load Makie as well (see
[`nomogramplot`](@ref)).

# Keywords
- `isopleths`, `style`: forwarded to [`drawing`](@ref) when given a `Nomogram`.
- `pxperunit = 4.0`: pixels per chart unit; the only size control you normally need.
- `padding = 2.0`: margin in chart units, added outside the computed bounding box.
- `background = "white"`, `axiscolor`, `textcolor`, `isoplethcolors`: colours, as any
  string an SVG renderer accepts. Pass `background = "none"` for a transparent chart.
- `gridcolor`: the colour of a net's strands and their labels. By default they are drawn
  like the axes, since they are read like them; pass a lighter colour to make a dense net
  less prominent.
- `font`: CSS font stack for all text.
"""
function svg(
    d::NomogramDrawing;
    pxperunit::Real = 4.0,
    padding::Real = 2.0,
    background::AbstractString = "white",
    axiscolor::AbstractString = "#333333",
    gridcolor::Union{Nothing,AbstractString} = nothing,
    textcolor::AbstractString = "#1a1a1a",
    isoplethcolors = DEFAULT_ISOPLETH_COLORS,
    isoplethdashes = DEFAULT_ISOPLETH_DASHES,
    font::AbstractString = "Helvetica, Arial, system-ui, sans-serif",
)
    st = d.style
    xmin, ymin, xmax, ymax = d.bbox
    pad = float(padding)
    k = float(pxperunit)
    W = (xmax - xmin + 2pad) * k
    H = (ymax - ymin + 2pad) * k
    # Chart coordinates run y-up; SVG runs y-down.
    X(x) = (x - xmin + pad) * k
    Y(y) = (ymax + pad - y) * k

    io = IOBuffer()
    print(
        io,
        """<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="$(_num(W))" height="$(_num(H))" \
viewBox="0 0 $(_num(W)) $(_num(H))">
""",
    )
    if background != "none"
        println(io, """  <rect width="100%" height="100%" fill="$(_esc(background))"/>""")
    end
    println(io, """  <g font-family="$(_esc(font))" stroke-linecap="round">""")

    if !isempty(d.title)
        _text(
            io,
            "    ",
            X,
            Y,
            d.titlepos,
            d.title,
            st.titlesize * k,
            d.titlealign,
            0.0,
            textcolor,
            "600",
        )
    end

    for sd in d.scales
        println(io, """    <g id="scale-$(sd.name)">""")
        grid = sd.kind === :grid
        stroke = grid ? something(gridcolor, axiscolor) : axiscolor
        linewidth = (grid ? st.gridlinewidth : st.axislinewidth) * k
        pts = join((" $(_num(X(p[1]))),$(_num(Y(p[2])))" for p in sd.path), "")
        println(
            io,
            """      <polyline points="$(strip(pts))" fill="none" """ *
            """stroke="$(_esc(stroke))" stroke-width="$(_num(linewidth))"/>""",
        )
        for t in sd.ticks
            o = (
                t.position[1] + t.normal[1] * t.length,
                t.position[2] + t.normal[2] * t.length,
            )
            println(
                io,
                """      <line x1="$(_num(X(t.position[1])))" """ *
                """y1="$(_num(Y(t.position[2])))" x2="$(_num(X(o[1])))" """ *
                """y2="$(_num(Y(o[2])))" stroke="$(_esc(axiscolor))" """ *
                """stroke-width="$(_num(st.tickwidth * k))"/>""",
            )
            isnothing(t.label) && continue
            _text(
                io,
                "      ",
                X,
                Y,
                t.labelpos,
                t.label,
                st.ticklabelsize * k,
                t.align,
                t.rotation,
                textcolor,
            )
        end
        if !isempty(sd.label)
            _text(
                io,
                "      ",
                X,
                Y,
                sd.titlepos,
                sd.label,
                st.titlesize * k,
                sd.titlealign,
                0.0,
                grid ? something(gridcolor, textcolor) : textcolor,
            )
        end
        println(io, "    </g>")
    end

    # A halo is a stroke in the background colour, painted under the letters.
    halo = background == "none" ? "white" : background
    for t in d.texts
        extra =
            t.halo ?
            """ stroke="$(_esc(halo))" stroke-width="$(_num(0.6 * k))" """ *
            """stroke-linejoin="round" paint-order="stroke\"""" : ""
        _text(
            io,
            "    ",
            X,
            Y,
            t.pos,
            t.text,
            _textsize(t, st) * k,
            t.align,
            t.rotation,
            textcolor,
            "",
            extra,
        )
    end

    for (i, iso) in enumerate(d.isopleths)
        # A reading of a compound chart is several segments that belong together, so a
        # grouped isopleth takes its colour from the group rather than its position.
        c, dash =
            isopleth_style(iso.group == 0 ? i : iso.group, isoplethcolors, isoplethdashes)
        a, b = iso.segment
        println(io, """    <g id="isopleth-$i">""")
        println(
            io,
            """      <line x1="$(_num(X(a[1])))" y1="$(_num(Y(a[2])))" """ *
            """x2="$(_num(X(b[1])))" y2="$(_num(Y(b[2])))" stroke="$(_esc(c))" """ *
            """stroke-width="$(_num(st.isoplethwidth * k))"$(_dasharray(dash, k))/>""",
        )
        for p in iso.points
            println(
                io,
                """      <circle cx="$(_num(X(p[1])))" cy="$(_num(Y(p[2])))" """ *
                """r="$(_num(st.markersize * k))" fill="$(_esc(c))"/>""",
            )
        end
        if !isempty(iso.label)
            _text(
                io,
                "      ",
                X,
                Y,
                iso.labelpos,
                iso.label,
                st.ticklabelsize * k,
                iso.labelalign,
                0.0,
                c,
            )
        end
        println(io, "    </g>")
    end

    println(io, "  </g>")
    println(io, "</svg>")
    return String(take!(io))
end

function svg(
    nom::AbstractNomogram;
    style::NomogramStyle = NomogramStyle(),
    isopleths::AbstractVector{Isopleth} = Isopleth[],
    kwargs...,
)
    return svg(drawing(nom; style, isopleths); kwargs...)
end

"""
    save_svg(path, nom::AbstractNomogram; kwargs...) -> path
    save_svg(path, d::NomogramDrawing; kwargs...) -> path

Write [`svg`](@ref) output to a file. Keywords are those of `svg`.
"""
function save_svg(path::AbstractString, x; kwargs...)
    open(path, "w") do io
        write(io, svg(x; kwargs...))
    end
    return path
end

# Rich display: anything that can show an image -- a notebook, VS Code's plot pane, a
# Documenter @example block -- shows the chart itself rather than a description of it.
# Display a `drawing` to include isopleths: `drawing(nom; isopleths = [...])`.
Base.show(io::IO, ::MIME"image/svg+xml", x::Union{AbstractNomogram,NomogramDrawing}) =
    print(io, svg(x))
