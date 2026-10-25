```@meta
CurrentModule = Nomograms
```

# Bibliography

Nomography is a nineteenth-century subject with a long textbook tradition. These are the
works this documentation leans on, and where to go for more.

## The founding text

- Maurice d'Ocagne. *Traité de nomographie. Théorie des abaques, applications pratiques.*
  Paris: Gauthier-Villars, 1899. A public-domain scan is on the
  [Internet Archive](https://archive.org/details/traitdenomograp00ocaggoog).

  The alignment method — three points on one straight line, written as a vanishing
  determinant — is d'Ocagne's, and every chart in this package rests on it.

## Textbooks

- L. Ivan Epstein. *Nomography.* New York: Interscience, 1958.
- Edward Otto. *Nomography.* New York: Macmillan, 1963. Discusses the criterion of
  Saint-Robert and Massau's method for deciding when an equation can be charted.
- Douglas P. Adams. *Nomography: Theory and Application.* Hamden, CT: Archon, 1964.

## A modern overview

- Ron Doerfler. "On Jargon: The Lost Art of Nomography." *The UMAP Journal* 30(4), 2009,
  pp. 457–493. The author posts a [free PDF](https://deadreckonings.com/wp-content/uploads/2010/04/the_lost_art_of_nomography.pdf).

  The best single introduction. It derives the parallel-scale chart, the N-chart and the
  determinant method, works the curved-scale chart for the quadratic `w² + uw + v = 0`
  (after Adams and Otto), builds grid nomograms from determinants, and reshapes charts by
  projective transformation.

## When a formula has a chart

- P. de Saint-Robert. "De la résolution de certaines équations à trois variables par le
  moyen d'une règle glissante." *Memorie della Reale Accademia delle Scienze di Torino*
  XXV, 1871, pp. 53–72. Presented in 1867.

  The test that [`auto_nomogram`](@ref) applies — that ``z = F(x, y)`` can be put in the
  form ``h(z) = f(x) + g(y)`` exactly when
  ``\partial^2 \ln(F_x / F_y) / \partial x \, \partial y = 0`` — is attributed to
  Saint-Robert by V. A. Esaulov,
  ["The Art of Nomography"](https://kartsci.org/kocomu/computer-history/the-art-of-nomography/),
  and by Lucien Pirio,
  ["Géométrie anamorphique II"](https://web.archive.org/web/20230203100329/http://images.math.cnrs.fr/Geometrie-Anamorphique-II.html)
  in *Images des mathématiques* (CNRS), 16 April 2011, which notes it has been
  rediscovered many times since (linked as archived in 2023; the original address no
  longer responds reliably). We follow those two accounts; the 1871 paper itself was not
  consulted.

  The more general problem — which equations can be brought into a chartable form at all —
  was taken up by J. Massau, "Mémoire sur l'intégration graphique et ses applications",
  *Annales de l'Association des ingénieurs sortis des écoles spéciales de Gand* 7, 1884,
  pp. 53–132, and independently by Léon Lecornu, "Sur le problème de l'anamorphose",
  *Comptes rendus de l'Académie des sciences* 102, 1886, pp. 813–816.

## Related software

- [PyNomo](http://www.pynomo.org), by Leif Roschier, builds nomograms in Python from a
  description of each scale. P. Martínez-Pagán and L. Roschier collect PyNomo scripts for
  common engineering formulas in *Data in Brief*, 2022,
  [doi:10.1016/j.dib.2022.108661](https://doi.org/10.1016/j.dib.2022.108661). Nomograms.jl
  is an independent implementation and borrows no code from it.

## What is original here

The derivations in these pages were written for this package, starting from the
collinearity condition, rather than taken from the sources above; where a result is
classical, the sources above are where it comes from. Specific to this package:

- the test that a formula of three or more inputs splits into a chain — that it depends on
  a pair of inputs only through one combination of them exactly when ``F_x / F_y`` does not
  depend on the rest — and the recovery of every transform numerically from the formula
  itself;
- how a chain's stages are paired, mirrored and sized so that no axis lands on another;
- how ticks are chosen and thinned, and how labels are kept clear of each other and of the
  alignment lines;
- the isopleth palette, whose colour-difference figures come from the package's own check
  against colour-vision deficiency rather than from a published source.
