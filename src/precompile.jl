# Run the common paths once at precompile time, so a first chart in a fresh session draws
# in a second or two rather than spending ten compiling. Formulas reach the automatic
# constructions as `_Formula`, so the closures here stand in for any a user writes.
PrecompileTools.@setup_workload begin
    PrecompileTools.@compile_workload begin
        g = auto_nomogram(
            (r, n) -> (1 + r / 100)^n,
            :rate => (1, 20),
            :years => (1, 30),
            :growth,
        )
        evaluate(g; rate = 7, growth = 2)
        svg(drawing(g; isopleths = [isopleth(g; rate = 5, years = 10, label = "x")]))

        chain = auto_nomogram(
            (a, b, c) -> a * b / c,
            :a => (1, 10),
            :b => (1, 10),
            :c => (1, 10),
            :d,
        )
        svg(drawing(chain; isopleths = isopleth(chain; a = 2, b = 3, c = 4, label = "x")))
        # A chain one level deeper builds closures of its own; five inputs mostly reuses them.
        auto_nomogram(
            (a, b, c, d) -> a * b * c / d,
            :a => (1, 10),
            :b => (1, 10),
            :c => (1, 10),
            :d => (1, 10),
            :e,
        )

        nom = additive_nomogram(
            Variable(:a, ReciprocalScale(); domain = (1, 10)),
            Variable(:b, ReciprocalScale(); domain = (1, 10)),
            Variable(:c, ReciprocalScale()),
        )
        evaluate(nom; a = 4, b = 6)
        svg(nom; isopleths = [isopleth(nom; a = 4, b = 6)])
        svg(fan(nom))

        prod = multiplicative_nomogram(
            Variable(:a; domain = (1, 10)),
            Variable(:b; domain = (1, 10)),
            Variable(:c),
        )
        evaluate(prod; a = 2, b = 3)
        n = ratio_nomogram(
            Variable(:a; domain = (1, 10)),
            Variable(:b; domain = (1, 10)),
            Variable(:c),
        )
        svg(drawing(n; isopleths = [isopleth(n; a = 6, b = 3)]))

        quad = determinant_nomogram(
            Variable(:p; domain = (-6, 0)) => p -> (0.0, p),
            Variable(:q; domain = (0, 6)) => q -> (1.0, q),
            Variable(:z; domain = (0.2, 5)) => z -> (1 / (1 + z), -z^2 / (1 + z));
            relation = (p, q, z) -> z^2 + p * z + q,
        )
        svg(
            drawing(
                quad;
                isopleths = [isopleth(quad, evaluate_all(quad; p = -3, q = 2)[1])],
            ),
        )
        max_relation_residual(quad)

        conc = contour_chart(
            Variable(:x, LogScale(); domain = (1, 10)),
            Variable(:y; domain = (0, 1)),
            Variable(:z) => (x, y) -> x * (1 + y),
            Variable(:w) => (x, y) -> x - y;
            region = (x, y) -> y <= x,
            width = 40,
            height = 30,
        )
        evaluate(conc; z = 3, w = 2)
        svg(conc; crossings = [crossing(conc; x = 2, z = 3, label = "r")])
        svg(contour_chart((x, y) -> x * y, :x => (1, 10), :y => (1, 10)))
        disc = contour_chart(
            :u => (-1, 1),
            :v => (-1, 1),
            :r => (u, v) -> hypot(u, v),
            :s => (u, v) -> u + 2v;
            axes = false,
            region = (u, v) -> u^2 + v^2 <= 1,
            rings = [Variable(:a; domain = (-180, 180)) => t -> (1.1cosd(t), 1.1sind(t))],
            labelat = (r = :v => 0,),
            width = 30,
        )
        evaluate(disc; r = 0.5, a = 30)
        svg(disc; crossings = [crossing(disc, first(evaluate_all(disc; r = 0.5, s = 0.5)))])

        svg(
            project(
                nom,
                [(0, 0), (70, 0), (70, 100), (0, 100)] =>
                    [(14, 0), (56, 0), (70, 100), (0, 100)],
            ),
        )
    end
end
