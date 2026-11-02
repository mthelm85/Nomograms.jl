```@meta
CurrentModule = Nomograms
```

# Bibliography

Nomography dates from the nineteenth century and has a long textbook tradition. This page
lists the works this documentation draws on, and where to read more.

## The founding text

- Maurice d'Ocagne. *Traité de nomographie. Théorie des abaques, applications pratiques.*
  Paris: Gauthier-Villars, 1899. A public-domain scan is on the
  [Internet Archive](https://archive.org/details/traitdenomograp00ocaggoog).

  The alignment method, in which three points lie on one straight line and the condition
  is written as a vanishing determinant, is d'Ocagne's. Every chart in this package is
  based on it.

## Textbooks

- L. Ivan Epstein. *Nomography.* New York: Interscience, 1958.
- Edward Otto. *Nomography.* New York: Macmillan, 1963. Discusses the criterion of
  Saint-Robert and Massau's method for deciding when an equation can be charted.
- Douglas P. Adams. *Nomography: Theory and Application.* Hamden, CT: Archon, 1964.

## A modern overview

- Ron Doerfler. "On Jargon: The Lost Art of Nomography." *The UMAP Journal* 30(4), 2009,
  pp. 457–493. The author posts a [free PDF](https://deadreckonings.com/wp-content/uploads/2010/04/the_lost_art_of_nomography.pdf).

  This is the best single introduction. It derives the parallel-scale chart, the N-chart and the
  determinant method, works the curved-scale chart for the quadratic `w² + uw + v = 0`
  (after Adams and Otto), builds grid nomograms from determinants, and reshapes charts by
  projective transformation.

## When a formula has a chart

- P. de Saint-Robert. "De la résolution de certaines équations à trois variables par le
  moyen d'une règle glissante." *Memorie della Reale Accademia delle Scienze di Torino*
  XXV, 1871, pp. 53–72. Presented in 1867.

  [`auto_nomogram`](@ref) applies the test that ``z = F(x, y)`` can be put in the form
  ``h(z) = f(x) + g(y)`` exactly when
  ``\partial^2 \ln(F_x / F_y) / \partial x \, \partial y = 0``. Two accounts attribute
  this test to Saint-Robert: V. A. Esaulov,
  ["The Art of Nomography"](https://kartsci.org/kocomu/computer-history/the-art-of-nomography/),
  and Lucien Pirio,
  ["Géométrie anamorphique II"](https://web.archive.org/web/20230203100329/http://images.math.cnrs.fr/Geometrie-Anamorphique-II.html)
  in *Images des mathématiques* (CNRS), 16 April 2011. Pirio notes that the test has been
  rediscovered many times since. The link is to a copy archived in 2023, because the
  original address no longer responds reliably. We follow those two accounts; the 1871
  paper itself was not consulted.

  The more general problem, of which equations can be brought into a chartable form at
  all, was taken up by J. Massau, "Mémoire sur l'intégration graphique et ses
  applications", *Annales de l'Association des ingénieurs sortis des écoles spéciales de
  Gand* 7, 1884, pp. 53–132, and independently by Léon Lecornu, "Sur le problème de
  l'anamorphose", *Comptes rendus de l'Académie des sciences* 102, 1886, pp. 813–816.

## Charts read where curves cross

- Léon Lalanne. "Mémoire sur les tables graphiques et sur la géométrie anamorphique."
  *Annales des Ponts et Chaussées*, 1st semester 1846.

  The source of anamorphosis: choosing axis scales that straighten a chart's curves, which
  [`contour_chart`](@ref) does for a separable formula. Cited here from d'Ocagne's *Traité*,
  chapter II, which also gives the history of the names used on the
  [contour charts](35-contour-charts.md) page: Lalanne's *courbes d'égal élément*, the
  *courbes isoplèthes* of Vogler (*Sechs graphische Tafeln*, Berlin, 1877), and d'Ocagne's
  own *courbes cotées*. Lalanne's and Vogler's works were not consulted directly.

## The charts in the examples

These pages draw the following charts to show what the package can do. The charts and their
formulas belong to the works below.

- **Smith chart.** P. H. Smith. "Transmission Line Calculator." *Electronics* 12(1),
  January 1939, pp. 29–31. The chart is reported to have been proposed independently by
  Tōsaku Mizuhashi in Japan in 1937.
- **Moody diagram.** L. F. Moody. "Friction Factors for Pipe Flow." *Transactions of the
  ASME* 66(8), 1944, pp. 671–684. It plots the equation of C. F. Colebrook, "Turbulent Flow
  in Pipes, with Particular Reference to the Transition Region between the Smooth and Rough
  Pipe Laws", *Journal of the Institution of Civil Engineers* 11, 1939, pp. 133–156.
- **Psychrometric chart.** The saturation pressure uses the coefficients of O. A. Alduchov
  and R. E. Eskridge, "Improved Magnus Form Approximation of Saturation Vapor Pressure",
  *Journal of Applied Meteorology* 35(4), 1996, pp. 601–609.
- **Wulff net.** After G. Wulff, "Untersuchungen im Gebiete der optischen Eigenschaften
  isomorpher Krystalle", *Zeitschrift für Krystallographie und Mineralogie* 36, 1902.

## Related software

- [PyNomo](http://www.pynomo.org), by Leif Roschier, builds nomograms in Python from a
  description of each scale. P. Martínez-Pagán and L. Roschier collect PyNomo scripts for
  common engineering formulas in *Data in Brief*, 2022,
  [doi:10.1016/j.dib.2022.108661](https://doi.org/10.1016/j.dib.2022.108661). Nomograms.jl
  is an independent implementation and borrows no code from it.

## What is original here

The derivations in these pages were written for this package, starting from the
collinearity condition, and are not copied from the sources above. Where a result is
classical, it comes from those sources. The following are specific to this package:

- the test that a formula of three or more inputs splits into a chain (it depends on a
  pair of inputs only through one combination of them exactly when ``F_x / F_y`` does not
  depend on the rest), and the numerical recovery of every transform from the formula
  itself;
- how a chain's stages are paired, mirrored and sized so that no axis lands on another;
- the search for a chart with one curved scale, which fits the condition on `ln|G_u / G_v|`
  numerically along level sets of the curved variable (the form it looks for, with two
  scales straight and one curved, is classical);
- how a contour chart's curves are thinned, from the spacing the field's gradient gives at
  each point;
- how ticks are chosen and thinned, and how labels are kept clear of each other and of the
  alignment lines;
- the isopleth palette, whose colour-difference figures come from the package's own check
  against colour-vision deficiency rather than from a published source.
