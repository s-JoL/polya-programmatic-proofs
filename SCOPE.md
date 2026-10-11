# Scope of verification

`FullTorus.lean` proves polynomial identities, three strictly positive Bernstein expansions, positivity for every (x,y)∈[-1,1]² and hence every pair of real phases. Its smallest coefficients are 240, 313600 and 240844800. It does not encode raw finite-element integration, complex matrix determinant identification, Sylvester's criterion, the Fourier assembly or the full matrix inverse/inertia argument. Those steps are handwritten in the short proof.

`K3Interfaces.lean` proves 12π<377/10, exact rational parameter identities, selection between small/large strip mass, coverage by four analytic gates and a scalar Schur square identity. The analytic gates are theorem hypotheses. It does not define planar domains or Dirichlet eigenvalues and does not prove the unconditional PDE theorem.

The exact/interval certificate programs validate the numerical inputs to the k=3 argument. Their mathematical sufficiency depends on the analytic certificate criteria and reductions in `paper/k3.md`. Replay, negative controls and internal reconstruction are scoped verification, not external peer review or a formal proof of all of spectral analysis.

This release does not supply a Pólya counterexample, a proof for every k, a forbidden-node CJ bound or an actual-area return for arbitrary geometry.
