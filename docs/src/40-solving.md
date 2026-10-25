```@meta
CurrentModule = Nomograms
```

# Solving

You can read a chart numerically as well as by eye. [`evaluate`](@ref) does what you would
do with a straightedge: it finds the line through the known points, and where that line
crosses the unknown scale.

```@example solve
using Nomograms

nom = additive_nomogram(
    Variable(:R1, ReciprocalScale(); domain = (1, 10)),
    Variable(:R2, ReciprocalScale(); domain = (1, 10)),
    Variable(:R, ReciprocalScale()),
)
evaluate(nom; R1 = 4, R2 = 6)
```

Any variable can be the unknown, so you solve the inverse problem the same way:

```@example solve
evaluate(nom; R1 = 4, R = 2.4)
```

[`isopleth`](@ref) takes the same arguments and returns the line itself, ready to draw. To
draw a reading, pass it a solution from `evaluate`:

```@example solve
drawing(nom; isopleths = [isopleth(nom, evaluate(nom; R1 = 4, R = 2.4); label = "R₂ = 6")])
```

## More than one crossing

A straight line can meet a *curved* scale more than once, and each crossing is a valid
solution of the relation. The usual shortcut is to check the sign of the residual at the
two ends of the domain and bisect. That fails in two ways. With an even number of
crossings, both ends lie on the same side of the line, so it reports no solution. With an
odd number greater than one, it returns one of them arbitrarily and does not tell you
there are others.

Instead, [`solve_all_on`](@ref) searches the whole domain. [`evaluate`](@ref) still returns
a single answer. When there are several, it throws an error that lists them instead of
picking one. [`evaluate_all`](@ref) returns every solution:

```@example solve
quad = determinant_nomogram(
    Variable(:p; domain = (-6, 0)) => p -> (0.0, p),
    Variable(:q; domain = (0, 6)) => q -> (1.0, q),
    Variable(:z; domain = (0.2, 5)) => z -> (1 / (1 + z), -z^2 / (1 + z));
    relation = (p, q, z) -> z^2 + p * z + q,
)
try
    evaluate(quad; p = -3, q = 2)
catch err
    println(sprint(showerror, err))
end
```

```@example solve
[sol.z for sol in evaluate_all(quad; p = -3, q = 2)]
```

On a straight scale the residual is an affine function of position, and so monotone,
which means there is at most one root. Charts made only of straight scales never meet
this case. Searching the whole domain costs around 120 µs per `evaluate` call, against a
few microseconds for bracketing, which makes no practical difference when reading charts.

## Checking a chart

[`max_relation_residual`](@ref) solves the chart over a grid of inputs and compares every
answer with the chart's `relation`. For a chart built by [`auto_nomogram`](@ref), the
relation is the function you passed in. For other charts, it is the separable form the
chart was built from, or whatever you passed as `relation`.

```@example solve
max_relation_residual(quad)
```

## Units

If you load [Unitful](https://github.com/PainterQubits/Unitful.jl), a variable can have a
unit, taken from its domain or given with the `unit` keyword. You can then pass values in
any unit of the right dimension. Answers come back in the axis's own unit, and the unit is
printed after the axis title:

```@example units
using Nomograms, Unitful

room = multiplicative_nomogram(
    Variable(:w; domain = (1u"ft", 20u"ft"), label = "width"),
    Variable(:l; domain = (1u"ft", 20u"ft"), label = "length"),
    Variable(:A; unit = u"ft^2", label = "area");
    title = "Floor area",
)
evaluate(room; w = 3u"m", l = 4u"m")
```

```@example units
drawing(room; isopleths = [isopleth(room; w = 10u"ft", l = 12u"ft", label = "120 ft²")])
```

A plain number is taken to be in the axis's unit, so `evaluate(room; w = 10, l = 12)`
means feet. A value with the wrong dimension is refused with an error that names the
variable, and so is a quantity given to a variable that has no unit. Temperatures are
converted with their offsets, so a °C axis reads 68 °F as 20 °C.

Units are most useful with [`auto_nomogram`](@ref). Your formula is called with actual
quantities, so it can compute in any units it likes. Its output is converted to the
result's unit, which is taken from the formula itself if you do not give one. That
conversion checks the formula's dimensions before any chart is built:

```@example units
try
    auto_nomogram(
        (m, h) -> m / h,   # should be m / h^2
        Variable(:mass; domain = (40u"kg", 150u"kg")),
        Variable(:height; domain = (1.4u"m", 2.1u"m")),
        Variable(:BMI; unit = u"kg/m^2"),
    )
catch err
    println(sprint(showerror, err))
end
```

Internally, the chart stores plain numbers in each variable's unit. Ticks are chosen in
that unit, and the `values` of an [`Isopleth`](@ref) are in it too. Formulas involving
temperature need kelvin, because Unitful does not multiply °C or °F.
