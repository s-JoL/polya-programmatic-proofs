# Pólya-related programmatic proofs

Research preprint and reproducible proof companion by **s-JoL**, 2026-10-11.

This release contains two results at different verification levels:

| Result | Mathematical scope | Machine-checked scope |
|---|---|---|
| Twisted tensor Q2 torus / complete Q1 inverse compression | All integer sizes L≥1, real Bloch phases and 0<E≤8, assuming the full pencil is invertible | Exact three-minor factorizations, Bernstein identities and positivity on the entire cosine square, hence all real phases |
| Third planar Dirichlet eigenvalue | Manuscript argument: &#124;Ω&#124;λ₃(Ω)>12π for every planar open set of finite positive area, multiplicities included | Exact/Arb B,T,S,N certificates and constants; Lean checks explicit rational interfaces, branch coverage and conditional selection logic |

The complete PDE theorem and complete finite-matrix inertia return are **not fully formalized in Lean**. The analytic arguments are in the manuscripts. The all-domain, all-k Pólya problem remains unsolved. No claim of priority or journal acceptance is made. AI assisted derivation, drafting and reconstruction; passing code alone does not establish the unformalized mathematical return.

## Read the proofs

- [Preprint overview source](paper/preprint.tex) (PDF compilation not yet confirmed)
- [Full torus proof](paper/torus.md)
- [Analytic certificate criteria](paper/certificate-criteria.md)
- [Short proofs in Chinese](paper/short-proofs.zh.md)
- [Complete k=3 analytic manuscript](paper/k3.md)
- [Classical assumptions and references](paper/classical-inputs.md)
- [Precise Lean coverage](SCOPE.md)

## Reproduce

Use Python 3.12+ and the official [Lean/elan toolchain](https://lean-lang.org/install/).
The toolchain is fixed at Lean 4.34.0 and mathlib at commit `5ed2965256430c3649e86755f9576b54eca72435`.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
lake update
lake exe cache get
bash scripts/verify.sh
```

The k=3 replay takes several minutes; allow roughly 8 GiB address-space budget for Lean and sufficient disk space for mathlib. Verification writes only to ignored `build/` scratch space and never regenerates shipped certificate inputs. B and the torus checker use exact rational arithmetic; T/S/N use Arb interval arithmetic through python-flint, not a libm error assumption. All gates have negative controls, including a corrupted T node and a false Lean lower bound.

`lake build` builds the two public Lean modules. `scripts/check_lean.py` additionally recompiles each source and audits the printed axioms. Standard Lean axioms `propext`, `Classical.choice`, `Quot.sound` may occur; no custom mathematical axiom or sorry is used.

## References and status

Classical inequalities and rearrangement theory are credited in the manuscript; some textbook citations remain chapter-level and are marked as such. Related Pólya results include [Filonov–Levitin–Polterovich–Sher, Pólya's conjecture for Euclidean balls](https://arxiv.org/abs/2203.07696). Schur inertia is classical; see [Haynsworth (1968)](https://doi.org/10.1016/0024-3795(68)90050-5). This companion supplies exact, inspectable certificates rather than claiming these classical tools as new.

## License

The newly released companion code and manuscripts are provided under the MIT license. Lean and mathlib are external dependencies with their own licenses. See `THIRD_PARTY.md`.
