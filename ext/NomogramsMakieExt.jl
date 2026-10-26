module NomogramsMakieExt

using Makie: Makie, Figure, Axis, DataAspect, Point2f
using Makie: hidedecorations!, hidespines!, lines!, linesegments!, scatter!, text!
using Nomograms: Nomograms, AbstractNomogram, NomogramDrawing, NomogramStyle, Isopleth
using Nomograms: ConcurrencyChart, Reading
using Nomograms: drawing, nomogramplot, nomogramplot!

const DEFAULT_COLORS = Nomograms.DEFAULT_ISOPLETH_COLORS
const DEFAULT_DASHES = Nomograms.DEFAULT_ISOPLETH_DASHES

_p(p) = Point2f(p[1], p[2])

"""
    nomogramplot!(ax, d::NomogramDrawing; kwargs...)

Draw a laid-out nomogram into an existing Makie axis. Lengths in the drawing are data
units; `pxperunit` converts them to the pixel units Makie uses for text and strokes, and
should match the data-to-pixel ratio of `ax` for the chart to look as intended.
"""
function Nomograms.nomogramplot!(
    ax::Axis,
    d::NomogramDrawing;
    pxperunit::Real = 4.0,
    axiscolor = "#333333",
    gridcolor = nothing,
    textcolor = "#1a1a1a",
    isoplethcolors = DEFAULT_COLORS,
    isoplethdashes = DEFAULT_DASHES,
    font = Makie.automatic,
    titlefont = :bold,
)
    st = d.style
    k = float(pxperunit)

    if !isempty(d.title)
        text!(
            ax,
            _p(d.titlepos);
            text = d.title,
            fontsize = st.titlesize * k,
            align = d.titlealign,
            color = textcolor,
            font = titlefont,
        )
    end

    for sd in d.scales
        grid = sd.kind === :grid
        lines!(
            ax,
            [_p(p) for p in sd.path];
            color = grid ? something(gridcolor, axiscolor) : axiscolor,
            linewidth = (grid ? st.gridlinewidth : st.axislinewidth) * k,
        )
        segs = Point2f[]
        for t in sd.ticks
            push!(segs, _p(t.position))
            push!(
                segs,
                Point2f(
                    t.position[1] + t.normal[1] * t.length,
                    t.position[2] + t.normal[2] * t.length,
                ),
            )
        end
        isempty(segs) ||
            linesegments!(ax, segs; color = axiscolor, linewidth = st.tickwidth * k)
        for t in sd.ticks
            isnothing(t.label) && continue
            text!(
                ax,
                _p(t.labelpos);
                text = t.label,
                fontsize = st.ticklabelsize * k,
                align = t.align,
                rotation = t.rotation,
                color = textcolor,
                font = font,
            )
        end
        isempty(sd.label) && continue
        text!(
            ax,
            _p(sd.titlepos);
            text = sd.label,
            fontsize = st.titlesize * k,
            align = sd.titlealign,
            color = grid ? something(gridcolor, textcolor) : textcolor,
            font = font,
        )
    end

    for t in d.texts
        text!(
            ax,
            _p(t.pos);
            text = t.text,
            fontsize = st.titlesize * k,
            align = t.align,
            color = textcolor,
            font = font,
        )
    end

    for (i, iso) in enumerate(d.isopleths)
        c, dash = Nomograms.isopleth_style(
            iso.group == 0 ? i : iso.group,
            isoplethcolors,
            isoplethdashes,
        )
        a, b = iso.segment
        lines!(
            ax,
            [_p(a), _p(b)];
            color = c,
            linewidth = st.isoplethwidth * k,
            linestyle = dash,
        )
        scatter!(
            ax,
            [_p(p) for p in iso.points];
            color = c,
            markersize = 2 * st.markersize * k,
        )
        isempty(iso.label) && continue
        text!(
            ax,
            _p(iso.labelpos);
            text = iso.label,
            fontsize = st.ticklabelsize * k,
            align = iso.labelalign,
            color = c,
            font = font,
        )
    end
    return ax
end

function Nomograms.nomogramplot!(
    ax::Axis,
    nom::AbstractNomogram;
    style::NomogramStyle = NomogramStyle(),
    isopleths::AbstractVector{Isopleth} = Isopleth[],
    kwargs...,
)
    return Nomograms.nomogramplot!(ax, drawing(nom; style, isopleths); kwargs...)
end

"""
    nomogramplot(nom::AbstractNomogram; kwargs...) -> Figure

Draw a nomogram as a standalone Makie `Figure`, sized from the drawing's bounding box so
that one chart unit is `pxperunit` pixels. Save it with Makie's `save`, which gives PDF,
SVG or PNG depending on the extension.
"""
function Nomograms.nomogramplot(
    d::NomogramDrawing;
    pxperunit::Real = 4.0,
    padding::Real = 2.0,
    backgroundcolor = :white,
    kwargs...,
)
    xmin, ymin, xmax, ymax = d.bbox
    pad = float(padding)
    k = float(pxperunit)
    xlo, xhi = xmin - pad, xmax + pad
    ylo, yhi = ymin - pad, ymax + pad
    fig = Figure(;
        size = ((xhi - xlo) * k, (yhi - ylo) * k),
        figure_padding = 0,
        backgroundcolor,
    )
    ax = Axis(
        fig[1, 1];
        aspect = DataAspect(),
        limits = (xlo, xhi, ylo, yhi),
        backgroundcolor,
    )
    hidedecorations!(ax)
    hidespines!(ax)
    Nomograms.nomogramplot!(ax, d; pxperunit = k, kwargs...)
    return fig
end

function Nomograms.nomogramplot(
    nom::AbstractNomogram;
    style::NomogramStyle = NomogramStyle(),
    isopleths::AbstractVector{Isopleth} = Isopleth[],
    kwargs...,
)
    return Nomograms.nomogramplot(drawing(nom; style, isopleths); kwargs...)
end

# A concurrency chart is marked with readings, not isopleths.
function Nomograms.nomogramplot!(
    ax::Axis,
    c::ConcurrencyChart;
    style::NomogramStyle = NomogramStyle(),
    readings::AbstractVector{Reading} = Reading[],
    isopleths::AbstractVector{Isopleth} = Isopleth[],
    kwargs...,
)
    return Nomograms.nomogramplot!(ax, drawing(c; style, readings, isopleths); kwargs...)
end

function Nomograms.nomogramplot(
    c::ConcurrencyChart;
    style::NomogramStyle = NomogramStyle(),
    readings::AbstractVector{Reading} = Reading[],
    isopleths::AbstractVector{Isopleth} = Isopleth[],
    kwargs...,
)
    return Nomograms.nomogramplot(drawing(c; style, readings, isopleths); kwargs...)
end

end # module
