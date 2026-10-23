const DEFAULT_ISOPLETH_COLORS = ("#c1272d", "#1b6ca8", "#2e7d32", "#9c27b0", "#ef6c00")

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
        """fill="$(_esc(fill))"$bold$transform>$(_esc(text))</text>""",
    )
    return nothing
end

"""
    svg(nom::Nomogram; kwargs...) -> String
    svg(d::NomogramDrawing; kwargs...) -> String

Render a nomogram as a standalone SVG document, returned as a `String`.

This backend has no dependencies beyond the standard library, so the package draws
publication-quality vector output out of the box; load Makie as well if you want the chart
inside a `Figure` (see [`nomogramplot`](@ref)).

# Keywords
- `isopleths`, `style`: forwarded to [`drawing`](@ref) when given a `Nomogram`.
- `pxperunit = 4.0`: pixels per chart unit; the only size control you normally need.
- `padding = 2.0`: margin in chart units, added outside the computed bounding box.
- `background = "white"`, `axiscolor`, `textcolor`, `isoplethcolors`: colours, as any
  string an SVG renderer accepts. Pass `background = "none"` for a transparent chart.
- `font`: CSS font stack for all text.
"""
function svg(
    d::NomogramDrawing;
    pxperunit::Real = 4.0,
    padding::Real = 2.0,
    background::AbstractString = "white",
    axiscolor::AbstractString = "#1a1a1a",
    textcolor::AbstractString = "#1a1a1a",
    isoplethcolors = DEFAULT_ISOPLETH_COLORS,
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
        pts = join((" $(_num(X(p[1]))),$(_num(Y(p[2])))" for p in sd.path), "")
        println(
            io,
            """      <polyline points="$(strip(pts))" fill="none" """ *
            """stroke="$(_esc(axiscolor))" """ *
            """stroke-width="$(_num(st.axislinewidth * k))"/>""",
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
                textcolor,
                "600",
            )
        end
        println(io, "    </g>")
    end

    colors = collect(isoplethcolors)
    for (i, iso) in enumerate(d.isopleths)
        c = colors[mod1(i, length(colors))]
        a, b = iso.segment
        println(io, """    <g id="isopleth-$i">""")
        println(
            io,
            """      <line x1="$(_num(X(a[1])))" y1="$(_num(Y(a[2])))" """ *
            """x2="$(_num(X(b[1])))" y2="$(_num(Y(b[2])))" stroke="$(_esc(c))" """ *
            """stroke-width="$(_num(st.isoplethwidth * k))"/>""",
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
    nom::Nomogram;
    style::NomogramStyle = NomogramStyle(),
    isopleths::AbstractVector{Isopleth} = Isopleth[],
    kwargs...,
)
    return svg(drawing(nom; style, isopleths); kwargs...)
end

"""
    save_svg(path, nom::Nomogram; kwargs...) -> path
    save_svg(path, d::NomogramDrawing; kwargs...) -> path

Write [`svg`](@ref) output to a file. Keywords are those of `svg`.
"""
function save_svg(path::AbstractString, x; kwargs...)
    open(path, "w") do io
        write(io, svg(x; kwargs...))
    end
    return path
end
