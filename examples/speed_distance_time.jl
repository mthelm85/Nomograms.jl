# Speed, distance and time:
#
#     speed = distance / time
#
# A ratio, so this is an N-chart (or Z-chart): two vertical scales joined by a diagonal,
# with the result read off the diagonal.
#
# `multiplicative_nomogram` would also chart this, by taking logarithms. The N-chart earns
# its place by needing no logarithm: the two outer scales stay *linear*, so precision is
# even along them and a reader can interpolate between ticks by eye. That is what you want
# for everyday quantities.
#
# As with the parallel construction, `speed` declares no domain: the range it can reach is
# fixed by the other two and derived for you.

using Nomograms

nom = ratio_nomogram(
    Variable(:distance; domain = (1, 10), label = "distance (km)"),
    Variable(:time; domain = (1, 10), label = "time (h)"),
    Variable(:speed; label = "speed (km/h)");
    title = "speed = distance / time",
)

@show nom[:speed].domain                            # (0.1, 10.0)
@show evaluate(nom; distance = 6, time = 3).speed   # 2.0
@show evaluate(nom; distance = 9, speed = 2.25)     # ...and backwards: time = 4

# The diagonal is a projective scale, s = k*r / (k*r + 1), so it spreads a ratio across
# the page roughly logarithmically. Tick values are chosen to match: over a span of a
# decade or more the result axis takes decade-based ticks rather than linear ones, which
# would offer nothing at all below 2 — most of the diagonal. Pass `tickstyle = :linear`
# on the variable to override that.

isopleths = [
    isopleth(nom; distance = 6, time = 3, label = "6 / 3 = 2"),
    isopleth(nom; distance = 9, time = 4, label = "9 / 4 = 2.25"),
]

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)
save_svg(joinpath(outdir, "speed_distance_time.svg"), nom; isopleths)
