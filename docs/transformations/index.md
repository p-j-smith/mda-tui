# Transformations

There are currently five transformations available in MDAnalysis TUI:

- [translate all atoms by a vector](translate.md), which uses [`MDAnalysis.transformations.translate.translate`](https://docs.mdanalysis.org/stable/documentation_pages/transformations/translate.html#MDAnalysis.transformations.translate.translate)
- [center atoms or molecule in a box](center-in-box.md), which uses [`MDAnalysis.transformations.translate.center_in_box`](https://docs.mdanalysis.org/stable/documentation_pages/transformations/translate.html#MDAnalysis.transformations.translate.center_in_box)
- [wrap atoms in a box](wrap.md), which uses [`MDAnalysis.transformations.wrap.wrap`](https://docs.mdanalysis.org/stable/documentation_pages/transformations/wrap.html#MDAnalysis.transformations.wrap.wrap)
- [unwrap atom group fragments](unwrap.md), which uses [`MDAnalysis.transformations.wrap.unwrap`](https://docs.mdanalysis.org/stable/documentation_pages/transformations/wrap.html#MDAnalysis.transformations.wrap.unwrap)
- [prevent atoms from jumping across periodic boundaries](nojump.md), which uses [`MDAnalysis.transformations.nojump.NoJump`](https://docs.mdanalysis.org/stable/documentation_pages/transformations/nojump.html#MDAnalysis.transformations.nojump.NoJump)
