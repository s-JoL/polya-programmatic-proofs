import Mathlib.Analysis.Real.Pi.Bounds
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring

noncomputable section

namespace PolyaAudit

-- These are interface and coverage theorems. Analytic gates are explicit hypotheses.
-- The PDE, rearrangement and HJB certificates are not silently declared as axioms.
theorem spectral_threshold : 12 * Real.pi < (377 / 10 : ℝ) := by
  have h := Real.pi_lt_d4
  linarith

theorem rational_interfaces :
    2 * (9/20:ℝ) = 9/10 ∧
    1 - 4 * (1141/5000:ℝ) = 109/1250 ∧
    (109/1250:ℝ) < 879293/10000000 ∧
    (240/1000:ℝ) < 241/1000 ∧
    (240/1000:ℝ) * (9/20) = 108/1000 := by
  norm_num

theorem select_cut (q : ℝ) (balanced_good perpendicular_good : Prop)
    (small : q ≤ 2*(1141/5000:ℝ) → balanced_good)
    (large : 2*(1141/5000:ℝ) < q → perpendicular_good) :
    balanced_good ∨ perpendicular_good := by
  rcases le_or_gt q (2*(1141/5000:ℝ)) with h | h
  · exact Or.inl (small h)
  · exact Or.inr (large h)

theorem A_from_four_gates (lam mp w : ℝ)
    (_hm : 0 ≤ mp ∧ mp ≤ 1141/5000) (hw : 0 ≤ w)
    (neck : mp ≤ 3/100 → (458318/10000:ℝ) ≤ lam)
    (zero : 3/100 ≤ mp → w = 0 → (398/10:ℝ) ≤ lam)
    (thin : 3/100 ≤ mp → 0 < w → w ≤ 241/1000 → (377/10:ℝ) ≤ lam)
    (wide : 3/100 ≤ mp → 240/1000 ≤ w → (377/10:ℝ) < lam) :
    12 * Real.pi < lam := by
  have hp := spectral_threshold
  by_cases h0 : mp ≤ 3/100
  · have hn := neck h0
    linarith
  · have hm0 : 3/100 ≤ mp := le_of_lt (lt_of_not_ge h0)
    by_cases hw0 : w = 0
    · have hz := zero hm0 hw0
      linarith
    · have hwpos : 0 < w := lt_of_le_of_ne hw (Ne.symm hw0)
      by_cases ht : w ≤ 241/1000
      · have hh := thin hm0 hwpos ht
        exact lt_of_lt_of_le hp hh
      · have hs : 240/1000 ≤ w := by linarith
        have hh := wide hm0 hs
        exact lt_trans hp hh

-- The standard scalar Schur completion; no statement about geometry is encoded.
theorem schur_square (h k c u v : ℝ) (hc : c ≠ 0) :
    h*u^2 + 2*k*u*v + c*v^2 =
      (h-k^2/c)*u^2 + c*(v+k*u/c)^2 := by
  field_simp
  ring

#print axioms spectral_threshold
#print axioms rational_interfaces
#print axioms select_cut
#print axioms A_from_four_gates
#print axioms schur_square
end PolyaAudit


