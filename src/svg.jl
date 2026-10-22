const DEFAULT_ISOPLETH_COLORS = ("#c1272d", "#1b6ca8", "#2e7d32", "#9c27b0", "#ef6c00")

function _esc(s::AbstractString)
    out = replace(s, '&' => "&amp;", '<' => "&lt;", '>' => "&gt;")
    return replace(out, '"' => "&quot;")
end

_num(x::Real) = @sprintf("%.2f", x)

# Horizontal anchoring implied by a tick's outward normal.
function _anchor(n::Point)
    n[1] < -0.3 && return "end"
    n[1] > 0.3 && return "start"
    return "middle"
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
    println(
        io,
        """  <g font-family="$(_esc(font))" fill="$(_esc(textcolor))" """ *
        """stroke-linecap="round">""",
    )

    if !isempty(d.title)
        println(
            io,
            """    <text x="$(_num(X(d.titlepos[1])))" """ *
            """y="$(_num(Y(d.titlepos[2])))" font-size="$(_num(st.titlesize * k))" """ *
            """text-anchor="middle" font-weight="600">$(_esc(d.title))</text>""",
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
            ox = t.position[1] + t.normal[1] * t.length
            oy = t.position[2] + t.normal[2] * t.length
            println(
                io,
                """      <line x1="$(_num(X(t.position[1])))" """ *
                """y1="$(_num(Y(t.position[2])))" x2="$(_num(X(ox)))" """ *
                """y2="$(_num(Y(oy)))" stroke="$(_esc(axiscolor))" """ *
                """stroke-width="$(_num(st.tickwidth * k))"/>""",
            )
            isnothing(t.label) && continue
            lx = ox + t.normal[1] * st.ticklabelgap
            ly = oy + t.normal[2] * st.ticklabelgap
            println(
                io,
                """      <text x="$(_num(X(lx)))" y="$(_num(Y(ly)))" """ *
                """font-size="$(_num(st.ticklabelsize * k))" """ *
                """text-anchor="$(_anchor(t.normal))" """ *
                """dominant-baseline="central">$(_esc(t.label))</text>""",
            )
        end
        if !isempty(sd.label)
            println(
                io,
                """      <text x="$(_num(X(sd.titlepos[1])))" """ *
                """y="$(_num(Y(sd.titlepos[2])))" """ *
                """font-size="$(_num(st.titlesize * k))" text-anchor="middle" """ *
                """font-weight="600">$(_esc(sd.label))</text>""",
            )
        end
        println(io, "    </g>")
    end

    colors = collect(isoplethcolors)
    for (i, iso) in enumerate(d.isopleths)
        c = _esc(colors[mod1(i, length(colors))])
        a, b = iso.segment
        println(io, """    <g id="isopleth-$i">""")
        println(
            io,
            """      <line x1="$(_num(X(a[1])))" y1="$(_num(Y(a[2])))" """ *
            """x2="$(_num(X(b[1])))" y2="$(_num(Y(b[2])))" stroke="$c" """ *
            """stroke-width="$(_num(st.isoplethwidth * k))"/>""",
        )
        for p in iso.points
            println(
                io,
                """      <circle cx="$(_num(X(p[1])))" cy="$(_num(Y(p[2])))" """ *
                """r="$(_num(st.markersize * k))" fill="$c"/>""",
            )
        end
        if !isempty(iso.label)
            println(
                io,
                """      <text x="$(_num(X(b[1] + st.ticklabelgap)))" """ *
                """y="$(_num(Y(b[2])))" font-size="$(_num(st.ticklabelsize * k))" """ *
                """text-anchor="start" dominant-baseline="central" """ *
                """fill="$c">$(_esc(iso.label))</text>""",
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
