import Mathlib.Basic.Real.Basic
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.NormNum


set_option maxHeartbeats 4000000
set_option maxRecDepth 8000

noncomputable section

namespace PolyaAudit


def minor1 (x y : ℝ) : ℝ := 1550 + 415 * y^1 + 1225 * x^1 + 635 * x^1 * y^1 + 75 * x^2 + 60 * x^2 * y^1

def bernstein1 (x y : ℝ) : ℝ := 560 * (1 * x^0 * (1-x)^2 * 1 * y^0 * (1-y)^1) + 240 * (1 * x^0 * (1-x)^2 * 1 * y^1 * (1-y)^0) + 1120 * (2 * x^1 * (1-x)^1 * 1 * y^0 * (1-y)^1) + 1830 * (2 * x^1 * (1-x)^1 * 1 * y^1 * (1-y)^0) + 1740 * (1 * x^2 * (1-x)^0 * 1 * y^0 * (1-y)^1) + 3960 * (1 * x^2 * (1-x)^0 * 1 * y^1 * (1-y)^0)

theorem minor1_identity (x y : ℝ) : minor1 (2*x-1) (2*y-1) = bernstein1 x y := by unfold minor1 bernstein1; ring

theorem minor1_positive (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1) (hy : 0 ≤ y ∧ y ≤ 1) : 0 < minor1 (2*x-1) (2*y-1) := by

  have hxnn : 0 ≤ x := hx.1
  have hynn : 0 ≤ y := hy.1
  have hx0 : 0 ≤ 1-x := sub_nonneg.mpr hx.2

  have hy0 : 0 ≤ 1-y := sub_nonneg.mpr hy.2

  have hp : (1 * x^0 * (1-x)^2 * 1 * y^0 * (1-y)^1) + (1 * x^0 * (1-x)^2 * 1 * y^1 * (1-y)^0) + (2 * x^1 * (1-x)^1 * 1 * y^0 * (1-y)^1) + (2 * x^1 * (1-x)^1 * 1 * y^1 * (1-y)^0) + (1 * x^2 * (1-x)^0 * 1 * y^0 * (1-y)^1) + (1 * x^2 * (1-x)^0 * 1 * y^1 * (1-y)^0) = (1:ℝ) := by ring

  have hw0 : 0 ≤ (1 * x^0 * (1-x)^2 * 1 * y^0 * (1-y)^1) := by positivity

  have hw1 : 0 ≤ (1 * x^0 * (1-x)^2 * 1 * y^1 * (1-y)^0) := by positivity

  have hw2 : 0 ≤ (2 * x^1 * (1-x)^1 * 1 * y^0 * (1-y)^1) := by positivity

  have hw3 : 0 ≤ (2 * x^1 * (1-x)^1 * 1 * y^1 * (1-y)^0) := by positivity

  have hw4 : 0 ≤ (1 * x^2 * (1-x)^0 * 1 * y^0 * (1-y)^1) := by positivity

  have hw5 : 0 ≤ (1 * x^2 * (1-x)^0 * 1 * y^1 * (1-y)^0) := by positivity

  have hb : (240:ℝ) ≤ bernstein1 x y := by
    unfold bernstein1
    linarith only [hp, hw0, hw1, hw2, hw3, hw4, hw5]

  rw [minor1_identity]
  have hm : (0:ℝ) < 240 := by norm_num
  exact lt_of_lt_of_le hm hb

#print axioms minor1_positive

def minor2 (x y : ℝ) : ℝ := 2340000 + 2517000 * y^1 + 656500 * y^2 + 25500 * y^3 + 2517000 * x^1 + 3642600 * x^1 * y^1 + 1246900 * x^1 * y^2 + 66900 * x^1 * y^3 + 656500 * x^2 + 1246900 * x^2 * y^1 + 544600 * x^2 * y^2 + 42600 * x^2 * y^3 + 25500 * x^3 + 66900 * x^3 * y^1 + 42600 * x^3 * y^2 + 3600 * x^3 * y^3

def bernstein2 (x y : ℝ) : ℝ := 313600 * (1 * x^0 * (1-x)^3 * 1 * y^0 * (1-y)^3) + (1388800/3:ℝ) * (1 * x^0 * (1-x)^3 * 3 * y^1 * (1-y)^2) + 504000 * (1 * x^0 * (1-x)^3 * 3 * y^2 * (1-y)^1) + 417600 * (1 * x^0 * (1-x)^3 * 1 * y^3 * (1-y)^0) + (1388800/3:ℝ) * (3 * x^1 * (1-x)^2 * 1 * y^0 * (1-y)^3) + (2844800/3:ℝ) * (3 * x^1 * (1-x)^2 * 3 * y^1 * (1-y)^2) + (14392000/9:ℝ) * (3 * x^1 * (1-x)^2 * 3 * y^2 * (1-y)^1) + (7068800/3:ℝ) * (3 * x^1 * (1-x)^2 * 1 * y^3 * (1-y)^0) + 504000 * (3 * x^2 * (1-x)^1 * 1 * y^0 * (1-y)^3) + (14392000/9:ℝ) * (3 * x^2 * (1-x)^1 * 3 * y^1 * (1-y)^2) + (33344000/9:ℝ) * (3 * x^2 * (1-x)^1 * 3 * y^2 * (1-y)^1) + (21184000/3:ℝ) * (3 * x^2 * (1-x)^1 * 1 * y^3 * (1-y)^0) + 417600 * (1 * x^3 * (1-x)^0 * 1 * y^0 * (1-y)^3) + (7068800/3:ℝ) * (1 * x^3 * (1-x)^0 * 3 * y^1 * (1-y)^2) + (21184000/3:ℝ) * (1 * x^3 * (1-x)^0 * 3 * y^2 * (1-y)^1) + 15641600 * (1 * x^3 * (1-x)^0 * 1 * y^3 * (1-y)^0)

theorem minor2_identity (x y : ℝ) : minor2 (2*x-1) (2*y-1) = bernstein2 x y := by unfold minor2 bernstein2; ring

theorem minor2_positive (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1) (hy : 0 ≤ y ∧ y ≤ 1) : 0 < minor2 (2*x-1) (2*y-1) := by

  have hxnn : 0 ≤ x := hx.1
  have hynn : 0 ≤ y := hy.1
  have hx0 : 0 ≤ 1-x := sub_nonneg.mpr hx.2

  have hy0 : 0 ≤ 1-y := sub_nonneg.mpr hy.2

  have hp : (1 * x^0 * (1-x)^3 * 1 * y^0 * (1-y)^3) + (1 * x^0 * (1-x)^3 * 3 * y^1 * (1-y)^2) + (1 * x^0 * (1-x)^3 * 3 * y^2 * (1-y)^1) + (1 * x^0 * (1-x)^3 * 1 * y^3 * (1-y)^0) + (3 * x^1 * (1-x)^2 * 1 * y^0 * (1-y)^3) + (3 * x^1 * (1-x)^2 * 3 * y^1 * (1-y)^2) + (3 * x^1 * (1-x)^2 * 3 * y^2 * (1-y)^1) + (3 * x^1 * (1-x)^2 * 1 * y^3 * (1-y)^0) + (3 * x^2 * (1-x)^1 * 1 * y^0 * (1-y)^3) + (3 * x^2 * (1-x)^1 * 3 * y^1 * (1-y)^2) + (3 * x^2 * (1-x)^1 * 3 * y^2 * (1-y)^1) + (3 * x^2 * (1-x)^1 * 1 * y^3 * (1-y)^0) + (1 * x^3 * (1-x)^0 * 1 * y^0 * (1-y)^3) + (1 * x^3 * (1-x)^0 * 3 * y^1 * (1-y)^2) + (1 * x^3 * (1-x)^0 * 3 * y^2 * (1-y)^1) + (1 * x^3 * (1-x)^0 * 1 * y^3 * (1-y)^0) = (1:ℝ) := by ring

  have hw0 : 0 ≤ (1 * x^0 * (1-x)^3 * 1 * y^0 * (1-y)^3) := by positivity

  have hw1 : 0 ≤ (1 * x^0 * (1-x)^3 * 3 * y^1 * (1-y)^2) := by positivity

  have hw2 : 0 ≤ (1 * x^0 * (1-x)^3 * 3 * y^2 * (1-y)^1) := by positivity

  have hw3 : 0 ≤ (1 * x^0 * (1-x)^3 * 1 * y^3 * (1-y)^0) := by positivity

  have hw4 : 0 ≤ (3 * x^1 * (1-x)^2 * 1 * y^0 * (1-y)^3) := by positivity

  have hw5 : 0 ≤ (3 * x^1 * (1-x)^2 * 3 * y^1 * (1-y)^2) := by positivity

  have hw6 : 0 ≤ (3 * x^1 * (1-x)^2 * 3 * y^2 * (1-y)^1) := by positivity

  have hw7 : 0 ≤ (3 * x^1 * (1-x)^2 * 1 * y^3 * (1-y)^0) := by positivity

  have hw8 : 0 ≤ (3 * x^2 * (1-x)^1 * 1 * y^0 * (1-y)^3) := by positivity

  have hw9 : 0 ≤ (3 * x^2 * (1-x)^1 * 3 * y^1 * (1-y)^2) := by positivity

  have hw10 : 0 ≤ (3 * x^2 * (1-x)^1 * 3 * y^2 * (1-y)^1) := by positivity

  have hw11 : 0 ≤ (3 * x^2 * (1-x)^1 * 1 * y^3 * (1-y)^0) := by positivity

  have hw12 : 0 ≤ (1 * x^3 * (1-x)^0 * 1 * y^0 * (1-y)^3) := by positivity

  have hw13 : 0 ≤ (1 * x^3 * (1-x)^0 * 3 * y^1 * (1-y)^2) := by positivity

  have hw14 : 0 ≤ (1 * x^3 * (1-x)^0 * 3 * y^2 * (1-y)^1) := by positivity

  have hw15 : 0 ≤ (1 * x^3 * (1-x)^0 * 1 * y^3 * (1-y)^0) := by positivity

  have hb : (313600:ℝ) ≤ bernstein2 x y := by
    unfold bernstein2
    linarith only [hp, hw0, hw1, hw2, hw3, hw4, hw5, hw6, hw7, hw8, hw9, hw10, hw11, hw12, hw13, hw14, hw15]

  rw [minor2_identity]
  have hm : (0:ℝ) < 313600 := by norm_num
  exact lt_of_lt_of_le hm hb

#print axioms minor2_positive

def minor3 (x y : ℝ) : ℝ := 5792715000 + 7699002750 * y^1 + 2915550000 * y^2 + 213250500 * y^3 + -46665000 * y^4 + -2477250 * y^5 + 7699002750 * x^1 + 11281713300 * x^1 * y^1 + 4582068900 * x^1 * y^2 + 384490200 * x^1 * y^3 + -74911650 * x^1 * y^4 + -5582700 * x^1 * y^5 + 2915550000 * x^2 + 4582068900 * x^2 * y^1 + 1624631200 * x^2 * y^2 + -64923400 * x^2 * y^3 + -85973200 * x^2 * y^4 + -5843100 * x^2 * y^5 + 213250500 * x^3 + 384490200 * x^3 * y^1 + -64923400 * x^3 * y^2 + -164601200 * x^3 * y^3 + -44295100 * x^3 * y^4 + -2749800 * x^3 * y^5 + -46665000 * x^4 + -74911650 * x^4 * y^1 + -85973200 * x^4 * y^2 + -44295100 * x^4 * y^3 + -9089800 * x^4 * y^4 + -559650 * x^4 * y^5 + -2477250 * x^5 + -5582700 * x^5 * y^1 + -5843100 * x^5 * y^2 + -2749800 * x^5 * y^3 + -559650 * x^5 * y^4 + -38700 * x^5 * y^5

def bernstein3 (x y : ℝ) : ℝ := 240844800 * (1 * x^0 * (1-x)^5 * 1 * y^0 * (1-y)^5) + 403128320 * (1 * x^0 * (1-x)^5 * 5 * y^1 * (1-y)^4) + 627630080 * (1 * x^0 * (1-x)^5 * 10 * y^2 * (1-y)^3) + 890071040 * (1 * x^0 * (1-x)^5 * 10 * y^3 * (1-y)^2) + 1104035840 * (1 * x^0 * (1-x)^5 * 5 * y^4 * (1-y)^1) + 1104691200 * (1 * x^0 * (1-x)^5 * 1 * y^5 * (1-y)^0) + 403128320 * (5 * x^1 * (1-x)^4 * 1 * y^0 * (1-y)^5) + 757075968 * (5 * x^1 * (1-x)^4 * 5 * y^1 * (1-y)^4) + 1315336192 * (5 * x^1 * (1-x)^4 * 10 * y^2 * (1-y)^3) + 2115106816 * (5 * x^1 * (1-x)^4 * 10 * y^3 * (1-y)^2) + 3129163776 * (5 * x^1 * (1-x)^4 * 5 * y^4 * (1-y)^1) + 4248289280 * (5 * x^1 * (1-x)^4 * 1 * y^5 * (1-y)^0) + 627630080 * (10 * x^2 * (1-x)^3 * 1 * y^0 * (1-y)^5) + 1315336192 * (10 * x^2 * (1-x)^3 * 5 * y^1 * (1-y)^4) + 2474830848 * (10 * x^2 * (1-x)^3 * 10 * y^2 * (1-y)^3) + 4266575104 * (10 * x^2 * (1-x)^3 * 10 * y^3 * (1-y)^2) + 6785337344 * (10 * x^2 * (1-x)^3 * 5 * y^4 * (1-y)^1) + 10034104320 * (10 * x^2 * (1-x)^3 * 1 * y^5 * (1-y)^0) + 890071040 * (10 * x^3 * (1-x)^2 * 1 * y^0 * (1-y)^5) + 2115106816 * (10 * x^3 * (1-x)^2 * 5 * y^1 * (1-y)^4) + 4266575104 * (10 * x^3 * (1-x)^2 * 10 * y^2 * (1-y)^3) + 7709332992 * (10 * x^3 * (1-x)^2 * 10 * y^3 * (1-y)^2) + 12719142912 * (10 * x^3 * (1-x)^2 * 5 * y^4 * (1-y)^1) + 19417845760 * (10 * x^3 * (1-x)^2 * 1 * y^5 * (1-y)^0) + 1104035840 * (5 * x^4 * (1-x)^1 * 1 * y^0 * (1-y)^5) + 3129163776 * (5 * x^4 * (1-x)^1 * 5 * y^1 * (1-y)^4) + 6785337344 * (5 * x^4 * (1-x)^1 * 10 * y^2 * (1-y)^3) + 12719142912 * (5 * x^4 * (1-x)^1 * 10 * y^3 * (1-y)^2) + 21373714432 * (5 * x^4 * (1-x)^1 * 5 * y^4 * (1-y)^1) + 32794460160 * (5 * x^4 * (1-x)^1 * 1 * y^5 * (1-y)^0) + 1104691200 * (1 * x^5 * (1-x)^0 * 1 * y^0 * (1-y)^5) + 4248289280 * (1 * x^5 * (1-x)^0 * 5 * y^1 * (1-y)^4) + 10034104320 * (1 * x^5 * (1-x)^0 * 10 * y^2 * (1-y)^3) + 19417845760 * (1 * x^5 * (1-x)^0 * 10 * y^3 * (1-y)^2) + 32794460160 * (1 * x^5 * (1-x)^0 * 5 * y^4 * (1-y)^1) + 49446092800 * (1 * x^5 * (1-x)^0 * 1 * y^5 * (1-y)^0)

theorem minor3_identity (x y : ℝ) : minor3 (2*x-1) (2*y-1) = bernstein3 x y := by unfold minor3 bernstein3; ring

theorem minor3_positive (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1) (hy : 0 ≤ y ∧ y ≤ 1) : 0 < minor3 (2*x-1) (2*y-1) := by

  have hxnn : 0 ≤ x := hx.1
  have hynn : 0 ≤ y := hy.1
  have hx0 : 0 ≤ 1-x := sub_nonneg.mpr hx.2

  have hy0 : 0 ≤ 1-y := sub_nonneg.mpr hy.2

  have hp : (1 * x^0 * (1-x)^5 * 1 * y^0 * (1-y)^5) + (1 * x^0 * (1-x)^5 * 5 * y^1 * (1-y)^4) + (1 * x^0 * (1-x)^5 * 10 * y^2 * (1-y)^3) + (1 * x^0 * (1-x)^5 * 10 * y^3 * (1-y)^2) + (1 * x^0 * (1-x)^5 * 5 * y^4 * (1-y)^1) + (1 * x^0 * (1-x)^5 * 1 * y^5 * (1-y)^0) + (5 * x^1 * (1-x)^4 * 1 * y^0 * (1-y)^5) + (5 * x^1 * (1-x)^4 * 5 * y^1 * (1-y)^4) + (5 * x^1 * (1-x)^4 * 10 * y^2 * (1-y)^3) + (5 * x^1 * (1-x)^4 * 10 * y^3 * (1-y)^2) + (5 * x^1 * (1-x)^4 * 5 * y^4 * (1-y)^1) + (5 * x^1 * (1-x)^4 * 1 * y^5 * (1-y)^0) + (10 * x^2 * (1-x)^3 * 1 * y^0 * (1-y)^5) + (10 * x^2 * (1-x)^3 * 5 * y^1 * (1-y)^4) + (10 * x^2 * (1-x)^3 * 10 * y^2 * (1-y)^3) + (10 * x^2 * (1-x)^3 * 10 * y^3 * (1-y)^2) + (10 * x^2 * (1-x)^3 * 5 * y^4 * (1-y)^1) + (10 * x^2 * (1-x)^3 * 1 * y^5 * (1-y)^0) + (10 * x^3 * (1-x)^2 * 1 * y^0 * (1-y)^5) + (10 * x^3 * (1-x)^2 * 5 * y^1 * (1-y)^4) + (10 * x^3 * (1-x)^2 * 10 * y^2 * (1-y)^3) + (10 * x^3 * (1-x)^2 * 10 * y^3 * (1-y)^2) + (10 * x^3 * (1-x)^2 * 5 * y^4 * (1-y)^1) + (10 * x^3 * (1-x)^2 * 1 * y^5 * (1-y)^0) + (5 * x^4 * (1-x)^1 * 1 * y^0 * (1-y)^5) + (5 * x^4 * (1-x)^1 * 5 * y^1 * (1-y)^4) + (5 * x^4 * (1-x)^1 * 10 * y^2 * (1-y)^3) + (5 * x^4 * (1-x)^1 * 10 * y^3 * (1-y)^2) + (5 * x^4 * (1-x)^1 * 5 * y^4 * (1-y)^1) + (5 * x^4 * (1-x)^1 * 1 * y^5 * (1-y)^0) + (1 * x^5 * (1-x)^0 * 1 * y^0 * (1-y)^5) + (1 * x^5 * (1-x)^0 * 5 * y^1 * (1-y)^4) + (1 * x^5 * (1-x)^0 * 10 * y^2 * (1-y)^3) + (1 * x^5 * (1-x)^0 * 10 * y^3 * (1-y)^2) + (1 * x^5 * (1-x)^0 * 5 * y^4 * (1-y)^1) + (1 * x^5 * (1-x)^0 * 1 * y^5 * (1-y)^0) = (1:ℝ) := by ring

  have hw0 : 0 ≤ (1 * x^0 * (1-x)^5 * 1 * y^0 * (1-y)^5) := by positivity

  have hw1 : 0 ≤ (1 * x^0 * (1-x)^5 * 5 * y^1 * (1-y)^4) := by positivity

  have hw2 : 0 ≤ (1 * x^0 * (1-x)^5 * 10 * y^2 * (1-y)^3) := by positivity

  have hw3 : 0 ≤ (1 * x^0 * (1-x)^5 * 10 * y^3 * (1-y)^2) := by positivity

  have hw4 : 0 ≤ (1 * x^0 * (1-x)^5 * 5 * y^4 * (1-y)^1) := by positivity

  have hw5 : 0 ≤ (1 * x^0 * (1-x)^5 * 1 * y^5 * (1-y)^0) := by positivity

  have hw6 : 0 ≤ (5 * x^1 * (1-x)^4 * 1 * y^0 * (1-y)^5) := by positivity

  have hw7 : 0 ≤ (5 * x^1 * (1-x)^4 * 5 * y^1 * (1-y)^4) := by positivity

  have hw8 : 0 ≤ (5 * x^1 * (1-x)^4 * 10 * y^2 * (1-y)^3) := by positivity

  have hw9 : 0 ≤ (5 * x^1 * (1-x)^4 * 10 * y^3 * (1-y)^2) := by positivity

  have hw10 : 0 ≤ (5 * x^1 * (1-x)^4 * 5 * y^4 * (1-y)^1) := by positivity

  have hw11 : 0 ≤ (5 * x^1 * (1-x)^4 * 1 * y^5 * (1-y)^0) := by positivity

  have hw12 : 0 ≤ (10 * x^2 * (1-x)^3 * 1 * y^0 * (1-y)^5) := by positivity

  have hw13 : 0 ≤ (10 * x^2 * (1-x)^3 * 5 * y^1 * (1-y)^4) := by positivity

  have hw14 : 0 ≤ (10 * x^2 * (1-x)^3 * 10 * y^2 * (1-y)^3) := by positivity

  have hw15 : 0 ≤ (10 * x^2 * (1-x)^3 * 10 * y^3 * (1-y)^2) := by positivity

  have hw16 : 0 ≤ (10 * x^2 * (1-x)^3 * 5 * y^4 * (1-y)^1) := by positivity

  have hw17 : 0 ≤ (10 * x^2 * (1-x)^3 * 1 * y^5 * (1-y)^0) := by positivity

  have hw18 : 0 ≤ (10 * x^3 * (1-x)^2 * 1 * y^0 * (1-y)^5) := by positivity

  have hw19 : 0 ≤ (10 * x^3 * (1-x)^2 * 5 * y^1 * (1-y)^4) := by positivity

  have hw20 : 0 ≤ (10 * x^3 * (1-x)^2 * 10 * y^2 * (1-y)^3) := by positivity

  have hw21 : 0 ≤ (10 * x^3 * (1-x)^2 * 10 * y^3 * (1-y)^2) := by positivity

  have hw22 : 0 ≤ (10 * x^3 * (1-x)^2 * 5 * y^4 * (1-y)^1) := by positivity

  have hw23 : 0 ≤ (10 * x^3 * (1-x)^2 * 1 * y^5 * (1-y)^0) := by positivity

  have hw24 : 0 ≤ (5 * x^4 * (1-x)^1 * 1 * y^0 * (1-y)^5) := by positivity

  have hw25 : 0 ≤ (5 * x^4 * (1-x)^1 * 5 * y^1 * (1-y)^4) := by positivity

  have hw26 : 0 ≤ (5 * x^4 * (1-x)^1 * 10 * y^2 * (1-y)^3) := by positivity

  have hw27 : 0 ≤ (5 * x^4 * (1-x)^1 * 10 * y^3 * (1-y)^2) := by positivity

  have hw28 : 0 ≤ (5 * x^4 * (1-x)^1 * 5 * y^4 * (1-y)^1) := by positivity

  have hw29 : 0 ≤ (5 * x^4 * (1-x)^1 * 1 * y^5 * (1-y)^0) := by positivity

  have hw30 : 0 ≤ (1 * x^5 * (1-x)^0 * 1 * y^0 * (1-y)^5) := by positivity

  have hw31 : 0 ≤ (1 * x^5 * (1-x)^0 * 5 * y^1 * (1-y)^4) := by positivity

  have hw32 : 0 ≤ (1 * x^5 * (1-x)^0 * 10 * y^2 * (1-y)^3) := by positivity

  have hw33 : 0 ≤ (1 * x^5 * (1-x)^0 * 10 * y^3 * (1-y)^2) := by positivity

  have hw34 : 0 ≤ (1 * x^5 * (1-x)^0 * 5 * y^4 * (1-y)^1) := by positivity

  have hw35 : 0 ≤ (1 * x^5 * (1-x)^0 * 1 * y^5 * (1-y)^0) := by positivity

  have hb : (240844800:ℝ) ≤ bernstein3 x y := by
    unfold bernstein3
    linarith only [hp, hw0, hw1, hw2, hw3, hw4, hw5, hw6, hw7, hw8, hw9, hw10, hw11, hw12, hw13, hw14, hw15, hw16, hw17, hw18, hw19, hw20, hw21, hw22, hw23, hw24, hw25, hw26, hw27, hw28, hw29, hw30, hw31, hw32, hw33, hw34, hw35]

  rw [minor3_identity]
  have hm : (0:ℝ) < 240844800 := by norm_num
  exact lt_of_lt_of_le hm hb

#print axioms minor3_positive



def p (x : ℝ) := 5*(2+x)
def q (x : ℝ) := (16-x-x^2)/2
def r (x : ℝ) := 6*(1-x)
def v (x : ℝ) := x^2+24*x+39
def d1 (x y : ℝ) := 5*(r y*q x+p y*v x)-8*p y*q x
def d2 (x y : ℝ) := 5*(v y*p x+q y*r x)-8*q y*p x
def d3 (x y : ℝ) := 5*(v y*q x+q y*v x)-8*q y*q x
def u (x y : ℝ) := (25/2:ℝ)*(3*(x+y)-10)
def f (x y : ℝ) := (25/2:ℝ)*v x+(15*y-35)*q x
def g (x y : ℝ) := (25/2:ℝ)*v y+(15*x-35)*q y

theorem minor1_factored (x y : ℝ) : minor1 x y = d1 x y := by
  unfold minor1 d1 p q r v; ring

theorem minor2_factored (x y : ℝ) :
    minor2 x y = d1 x y*d2 x y-4*(1+x)*(1+y)*(u x y)^2 := by
  unfold minor2 d1 d2 p q r v u; ring

theorem minor3_factored (x y : ℝ) : minor3 x y =
    d1 x y*d2 x y*d3 x y+8*(1+x)*(1+y)*u x y*f x y*g x y
    -2*(1+x)*d1 x y*(g x y)^2-2*(1+y)*d2 x y*(f x y)^2
    -4*(1+x)*(1+y)*d3 x y*(u x y)^2 := by
  simp only [minor3, d1, d2, d3, p, q, r, v, u, f, g]; ring

theorem all_cosine_box (x y : ℝ) (hx : -1 ≤ x ∧ x ≤ 1)
    (hy : -1 ≤ y ∧ y ≤ 1) :
    0 < minor1 x y ∧ 0 < minor2 x y ∧ 0 < minor3 x y := by
  have htx : 0 ≤ (x+1)/2 ∧ (x+1)/2 ≤ 1 := by constructor <;> linarith
  have hty : 0 ≤ (y+1)/2 ∧ (y+1)/2 ≤ 1 := by constructor <;> linarith
  have ex : 2*((x+1)/2)-1=x := by ring
  have ey : 2*((y+1)/2)-1=y := by ring
  have h1 := minor1_positive ((x+1)/2) ((y+1)/2) htx hty
  have h2 := minor2_positive ((x+1)/2) ((y+1)/2) htx hty
  have h3 := minor3_positive ((x+1)/2) ((y+1)/2) htx hty
  rw [ex,ey] at h1 h2 h3
  exact ⟨h1,h2,h3⟩

theorem all_real_phases (s t : ℝ) :
    0 < minor1 (Real.cos s) (Real.cos t) ∧
    0 < minor2 (Real.cos s) (Real.cos t) ∧
    0 < minor3 (Real.cos s) (Real.cos t) := by
  exact all_cosine_box _ _ ⟨Real.neg_one_le_cos s,Real.cos_le_one s⟩
    ⟨Real.neg_one_le_cos t,Real.cos_le_one t⟩

#print axioms minor1_factored
#print axioms minor2_factored
#print axioms minor3_factored
#print axioms all_real_phases

end PolyaAudit

