```@meta
CurrentModule = Nomograms
```

# Solving

A chart is also a calculator, and reading it numerically is the same operation as reading
it by eye: find the line through the known points, and where it crosses the unknown scale.

```@example solve
using Nomograms

nom = additive_nomogram(
    Variable(:R1, ReciprocalScale(); domain = (1, 10)),
    Variable(:R2, ReciprocalScale(); domain = (1, 10)),
    Variable(:R, ReciprocalScale()),
)
evaluate(nom; R1 = 4, R2 = 6)
```

Any variable may be the unknown, so the inverse problem needs nothing new:

```@example solve
evaluate(nom; R1 = 4, R = 2.4)
```

[`isopleth`](@ref) takes the same arguments and returns the line itself, ready to draw.
Pass it a solution from `evaluate` to draw that reading:

```@example solve
drawing(nom; isopleths = [isopleth(nom, evaluate(nom; R1 = 4, R = 2.4); label = "R₂ = 6")])
```

## More than one crossing

A straight line can meet a *curved* scale more than once, and each crossing is a genuine
solution of the relation. The usual shortcut — checking the sign of the residual at the two
ends of the domain and bisecting — gets this wrong twice over: it reports "no solution" for
an even number of crossings, because both ends lie on the same side of the line, and it
returns an arbitrary one of an odd number without saying so.

So [`solve_all_on`](@ref) searches the whole domain. [`evaluate`](@ref) still returns a
single answer, but throws rather than choosing when there are several, naming the values it
could not choose between; [`evaluate_all`](@ref) returns every solution:

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

A straight scale makes the residual an affine function of position, hence monotone, so
there is always at most one root: parallel-layout nomograms never encounter this. The cost
of searching rather than bracketing is around 120 µs per `evaluate` instead of a few, which
is immaterial for reading charts.

## Checking a chart

[`max_relation_residual`](@ref) solves on the geometry over a grid of inputs and compares
every answer with the chart's `relation`. For a chart built by [`auto_nomogram`](@ref) the
relation is the function you passed in; for the others it is the separable form the chart
was built from, or whatever you passed as `relation`.

```@example solve
max_relation_residual(quad)
```
