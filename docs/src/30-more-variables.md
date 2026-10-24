```@meta
CurrentModule = Nomograms
```

# Four or more variables

One alignment line can relate only three scales. Past that there are two ways on: chain
charts together through a shared scale, or let one point carry two variables at once.
[`auto_nomogram`](@ref) builds chains by itself for formulas that allow them; this page is
about building either by hand.

## Chains and turning lines

A fourth variable means two stages sharing an intermediate, and [`compound`](@ref)
positions them so the shared scale coincides.

```@example more
using Nomograms

stage1 = multiplicative_nomogram(
    Variable(:L; domain = (1, 10), label = "length"),
    Variable(:W; domain = (1, 10), label = "width"),
    Variable(:t; label = "", ticks = Float64[]);
    width = 40,
)
stage2 = multiplicative_nomogram(
    Variable(:t; domain = stage1[:t].domain, label = "", ticks = Float64[]),
    Variable(:H; domain = (1, 10), label = "height"),
    Variable(:V; label = "volume");
    width = 130,
)
box = compound(stage1, stage2; title = "V = L × W × H")
drawing(box; isopleths = isopleth(box; L = 3, W = 4, H = 2, label = "3 × 4 × 2"))
```

The reader lays a line across the first stage, notes where it crosses the shared scale,
then pivots about that point and reads the second — which is why that scale is called the
turning line. Giving it an empty tick list draws it bare, as here; leave the ticks on and
it carries the intermediate value, which is worth doing when the intermediate means
something.

Solving propagates rather than following a fixed direction: any stage with exactly one
unknown is read, its answer becomes known, and the pass repeats. So inputs may be given in
whatever combination closes the chain,

```@example more
evaluate(box; L = 3, W = 4, V = 24).H
```

and a combination that cannot close it is refused, naming what stayed unknown. Supplying
`t`, `H` and `V`, for instance, pins the product `L × W` but neither factor.

Alignment is a similarity — rotate, scale, translate — carrying the later stage's shared
scale onto the placed one, checked along the whole domain rather than only at the two
points that fixed it. Being affine, it cannot change what either stage computes.

!!! note "The shared scale must be the same scale in both stages"
    Both stages above are multiplicative, so both place `t` logarithmically and the two
    coincide. Splitting `V = (L + W) × H` the same way would not work: the first stage
    puts `t` on a linear scale and the second on a logarithmic one, and no amount of
    moving the stages reconciles that. `compound` checks and says so rather than drawing
    a chart that lies.

Stages interleave rather than sitting side by side, because a parallel-axis construction
puts its result *between* its inputs: one input of the first stage always lands inside the
second. Sizing the stages decides where — here the first is 40 wide and the second 130, so
`width` falls between the turning line and `volume` rather than on top of `volume`.
`auto_nomogram` makes that choice itself, and [`project`](@ref) will reshape the result if
the arrangement is still awkward.

## When a pivot will not do: nets

A pivot works only when the first two variables collapse into a single intermediate. When
they do not, one of the three points an alignment line crosses must carry two variables at
once. That is a [`Net`](@ref): a two-parameter surface `(u, v) -> (x, y)`, drawn as a mesh —
the grid nomogram of the classic texts, which goes back to d'Ocagne and which Doerfler builds
from determinants (see the [Bibliography](96-references.md)).

Deriving one follows the same pattern as a curved scale. With `a` on `x = 0` at height
`f(a)` and `b` on `x = 1` at height `g(b)`, a net point is crossed exactly when

```math
\alpha f(a) + \beta g(b) = \gamma, \qquad
x_3 = \frac{\beta}{\alpha + \beta}, \quad y_3 = \frac{\gamma}{\alpha + \beta}
```

so any relation whose two net variables enter as the coefficients has a net. The general
quadratic `a·x² + b·x + c = 0`, written `a·x² + b·x = -c`, has `α = x²`, `β = x`, `γ = -c`:

```@example more
gq = net_nomogram(
    Net(
        Variable(:x; domain = (0.5, 5)),
        Variable(:c; domain = (0, 4)),
        (x, c) -> (1 / (1 + x), -c / (x * (1 + x))),
    ),
    Variable(:a; domain = (0.5, 3)) => a -> (0.0, a),
    Variable(:b; domain = (-6, 0)) => b -> (1.0, b);
    relation = (x, c, a, b) -> a * x^2 + b * x + c,
    title = "ax² + bx + c = 0",
)
roots = evaluate_all(gq; a = 1, b = -3, c = 2)
drawing(gq; isopleths = [isopleth(gq, roots[1]; label = "x² − 3x + 2 = 0")])
```

```@example more
[sol.x for sol in roots]
```

Which problem [`evaluate`](@ref) solves depends on what is unknown. With both net variables
known the net contributes a single point and it is the ordinary three-point alignment. With
one unknown, the line is already fixed by the two scales, and the answer is where it
crosses the strand for the known net variable — a curve, which the line may cross twice.
Here it does, because a quadratic has two roots.

The mesh is genuinely two-dimensional only when `(α : β : γ)` sweeps a region. When it
sweeps just a curve the net collapses onto a line, and the relation wanted a pivot instead —
which is a useful check on whether you have reached for the right construction.

Strands are scales of `kind = :grid`, drawn lighter than the scales a reader takes values
from. Each family is labelled at whichever end it fans out, since at the end where strands
converge the labels pile up and get nudged next to the wrong curve; set `titleend` on a
[`Variable`](@ref) to choose by hand.
