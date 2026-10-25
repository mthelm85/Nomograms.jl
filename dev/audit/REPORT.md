# Numerical audit

An audit of the numerical methods Nomograms.jl uses to build, check and solve charts,
run against the package at commit `73e0d88`. Items 1 and 2 were fixed after it. Every number below comes from a script in
this folder; `julia dev/audit/setup.jl` creates the environment, then run any `e*.jl`.

## Summary

| # | Area | Verdict | Action |
|---|---|---|---|
| 1 | Root finding on curved scales | **Bug**: missed some pairs of close crossings | **Fixed** |
| 2 | Text width estimate | **Bug**: could break the no-overlap guarantee | **Fixed** |
| 3 | Separability and pivot thresholds | Sound, with a wide margin | None |
| 4 | Finite differences vs automatic differentiation | FD adequate for today's features | Revisit with item 6 |
| 5 | The built-in accuracy check | Sound, but measures the wrong thing at the small end | Small improvement |
| 6 | Detecting curved-scale charts (prototype) | Feasible; separates cases cleanly | Design decision needed |

## 1. Root finding misses some pairs of close crossings — bug

`e5_roots.jl`, `e5b_miss.jl`, `e5c_fix.jl`. On the quadratic chart (roots on a curved scale
from 0.2 to 6), 20,000 random readings with two genuine roots:

- **7 were answered wrongly**: two roots 0.003 to 0.03 apart (0.05–0.5% of the scale's
  range, clearly separate on the chart) came back as one root or, in five cases, as **no
  solution**.
- Everything else is solid: roots 1e-6 apart are found separately to 3e-9; they merge only
  at 1e-8, below anything a chart can show. A double root counts once, lines just missing
  the curve find nothing, and roots at either end of the scale are found.

**Cause.** Between two close roots the residual (signed distance from the line) is
negative almost everywhere, with a narrow positive bump between the roots — in the first
failing case 0.013 wide on a 5.8-long scale, 2e-3 high against depths of −353.
`Roots.find_zeros` only finds a bump if one of its sample points lands in it. The 1e-9
padding `_search_bounds` adds shifts every sample, which is why the package misses cases
that a direct call happens to catch.

**More samples do not fix it**: raising `find_zeros`'s `no_pts` from 12 to 200 still misses
2 in 20,000, at 14× the cost.

**Recommended fix.** Sample the residual at points evenly spaced in *page distance* along
the curve, fine enough that two crossings any reader could tell apart cannot both fall
between samples, then refine each sign change with a bracketing solver; also check local
minima of |residual| for tangency. That turns "usually finds both" into a stated resolution
guarantee.

**Fixed.** `solve_all_on` now searches by subdivision with a proof of absence: the
residual is a cross product with the fixed vector `q − p`, so moving a page distance `d`
along the curve changes it by at most `|q − p|·d`, and a stretch whose end values exceed
that bound cannot hide a crossing. Stretches that might are halved down to 1e-7 of the
scale's length; a line that touches the curve there is a double root, located by
golden-section search. Afterwards: **0 misses in 20,000** readings (was 7); pairs 1e-4
apart are both found, accurate to 1e-11; closer than about 1e-5 of the scale's range they
count as one reading. A solve takes about 0.7 ms (was 0.4 ms).

## 2. The text width estimate can break the no-overlap guarantee — bug

`e6_text.jl`, `e6b_caps.jl`. Label boxes assume each character is `0.6 × font size` wide.
Measured against real font metrics:

| Text | Real / estimate (Arial = Helvetica metrics, as the SVG uses) |
|---|---|
| Digits, `1 3/8`, mixed-case titles | 0.74 – 0.93 — overestimated |
| `WAVELENGTH`, `MAXIMUM WATTAGE (W)` | 1.17, 1.05 — underestimated |
| `WWWW` | 1.57 |
| `iiii` | 0.37 |

Overestimating drops labels that would have fitted. Underestimating — any title heavy in
capitals or wide letters — lets two labels overlap that the layout believes are clear,
because the collision padding (0.3 units) is far smaller than the error.

**Recommended fix.** Embed a per-character advance table for Helvetica/Arial (ASCII plus the
symbols the package itself emits), keeping 0.6 em as the fallback for unknown characters. On
every string tested the table matched the real width exactly. No new dependency.

**Fixed.** `src/textwidth.jl` holds advance widths for 165 characters — printable ASCII
plus the symbols labels and units use — generated from the font by `gen_widths.jl`, with
0.6 em for anything else and a full em for East Asian wide characters.

## 3. Separability and pivot thresholds — sound

`e1_thresholds.jl`, `e1b_pivot.jl`.

- Separable formulas score ≤ 1e-9, including domains near 1e-7 and 1e7, six-decade ranges,
  `exp(10x)`, `x^10`, large offsets and scale factors. The 1e-6 threshold leaves a 1000×
  margin: **no wrong refusals**.
- Non-separable perturbations score in proportion to the cross term (≈ 9ε for
  `x + y + ε·x·y²`) and are accepted below ε ≈ 1e-7. When accepted, the chart's worst
  reading error is ≈ ε of the result's range (8e-8 at ε = 1e-7) — four orders of magnitude
  below reading precision.
- The chain (pivot) test behaves the same: combinable pairs ≤ 1e-9, others ≈ 1.

The threshold is, if anything, conservative: a formula refused at ε = 1e-6 would have made
a chart accurate to 1e-6. That headroom is what an "approximate chart" feature would use.

## 4. Finite differences vs automatic differentiation

`e2_ad.jl`.

- ForwardDiff lowers the separability test's noise floor from ~1e-10 to ~1e-16 and is
  slightly faster (2.0 µs vs 3.3 µs per test on a 9×9 grid).
- It does **not** improve the charts: recovered scales match the exact transform to ~3e-8
  either way; tabulation, not differentiation, sets that.
- It cannot run through `_Formula`, the `Float64`-only FunctionWrapper behind fast first
  calls. It does work through Unitful.

Not worth switching for today's features. It matters for item 6.

## 5. The built-in accuracy check

`e4_accuracy.jl`. Every auto chart is checked on a 9×9 (or 4ⁿ) grid, error relative to the
result's range. Against 3,000 random readings per chart:

- Dense sampling finds errors 1–20× larger than the grid reports — but every chart tested
  stays below 9e-5 of the axis length in page position: under 0.01 mm on a 10 cm axis,
  against roughly 0.1 mm for reading by eye.
- "Relative to range" understates error at the small end of a long-range result: compound
  growth scores 1e-9 by range but 7e-6 relative to the value.

**Recommended improvement** (low priority): measure the check in page position — what a
reader actually sees — and sample more points, randomly placed.

## 6. Detecting curved-scale charts — prototype

`e7_curved.jl`. A formula has a one-curved-scale chart with `w` on the curve exactly when
`L = ln|F_u/F_v|` splits as `φ(w) + a(u) − b(v)`. Along each level set `w = c` this means
`L_uv = A·w_u·w_v + B·w_uv` with constants `A = φ''(c)`, `B = φ'(c)`. The prototype fits `A`
and `B` on 21 level sets and checks the residual and `A = dB/dc`, for each choice of curved
variable:

| Formula | Passing assignment | Fit residual: pass / others | `A = dB/dc`: pass / others |
|---|---|---|---|
| `z = x·y` (separable) | all | 0 / 0 | ≤ 7e-8 |
| `z² + pz + q = 0` | `z` curved | 5e-16 / ≥ 0.14 | 2e-4 / ≥ 0.46 |
| `xz + y = z³` | `z` curved | 4e-16 / ≥ 0.02 | 2e-4 / ≥ 0.21 |
| `x² + xy + y = z` | **`x`** curved | 5e-16 / ≥ 0.03 | 6e-4 / ≥ 0.27 |
| `sin x + eʸ + xy = z` | none | — / ≥ 0.01 | — / ≥ 0.94 |

The test is feasible, separates the cases by many orders of magnitude, and finds
non-obvious assignments (the polynomial is linear in `y` and `z` with coefficients in `x`).

**Design decision needed.** The test needs exact third derivatives of an implicitly solved
variable, which automatic differentiation provides and finite differences would not do
reliably. Automatic differentiation cannot pass through the `Float64`-only `_Formula`
wrapper, so building this means either calling the user's formula unwrapped for detection
(slower first call) or a wrapper that also accepts dual numbers.

## Still worth asking the community

- The API conventions — argument order of `auto_nomogram`, plain numbers read in an axis's
  display unit, `drawing` vs `nomogramplot` — are about what Julia users expect, which no
  experiment settles.
- Literature on detecting nomographic forms (Massau, Lecornu, web geometry) that would
  improve on item 6.
