# 解析核对：k3.md 引用的经典结果——出处与假设核对

每条给出：在证明中的用处、标准出处、假设是否满足。证明中所有集合都是 ℝ² 的开集（或可测集），只假设有限测度，可以无界、不连通、边界任意不正则。
凡出处后标"（章节级）"的，我只核到章或节，没有逐字核对定理号，正式成文前应查原书。

---

## 1. Faber–Krahn 不等式

**内容。** 对任意开集 U ⊂ ℝ²，|U| < ∞：λ₁(U) ≥ κ/|U|，κ = πj₀₁²。
**H¹ 版本。** 若 v ∈ H¹(ℝ²) 在某个测度为 A 的可测集外 a.e. 为 0，则 ∫|∇v|² ≥ (κ/A)∫v²。

**出处。**
- G. Faber, Sitzungsber. Bayer. Akad. Wiss. München (1923) 169–172；E. Krahn, Math. Ann. 94 (1925) 97–100。
- A. Henrot, *Extremum Problems for Eigenvalues of Elliptic Operators*, Birkhäuser 2006, Thm 3.2.1。
- 无正则性、有限测度的版本：E. H. Lieb, M. Loss, *Analysis*, 2nd ed., AMS GSM 14 (2001), Lemma 7.17（对称递减重排不增动能）加 Ch. 3（重排保 Lᵖ 范数）。

**H¹ 版本的推导。** v* 支在面积 A 的闭圆盘 B̄ 内。由 Lieb–Loss Lemma 7.17，∫|∇v*|² ≤ ∫|∇v|²；该引理要求 v 在无穷远处消失，即 |{|v| > t}| < ∞，这里成立。圆盘是 Lipschitz 区域，所以 v* ∈ H¹₀(B)，再用 λ₁(B) = κ/A。

**用处与假设核对。**
| 位置 | 作用对象 | 假设 |
|---|---|---|
| 引理 2.2(iii)、引理 2.4 | 核 W_i、半块 W_i⁺：连续函数的正集，开，测度 ≤ 1/2 | ✓ 开、有限测度；不要求有界或正则 |
| 引理 4.2（非对称分支） | S_c 中落在半平面内的连通分支，开，面积 ≤ a | ✓ |
| 引理 4.4（帽质量） | u|_{C₊} 关于 y = τ 偶反射后的 H¹(ℝ²) 函数，在测度 ≤ 2α′ 的集合外为 0 | ✓ 用 H¹ 版本，不需要该集合开 |
| §4.6（w* = 0） | 迹为 0 的 v ∈ H¹₀(T) 零延拓，在 D*（测度 2m′）外为 0 | ✓ H¹ 版本 |

## 2. Krahn–Szegő 不等式

**内容。** λ₂(U) ≥ 2κ/|U|。

**出处。** E. Krahn, Acta Comm. Univ. Tartu (Dorpat) A9 (1926) 1–44；G. Pólya, Math. Z. 63 (1955) 331–337（其中给出 Szegő 的证明）；Henrot 2006, Thm 4.1.1。

**用处。** 证明中没有直接引用。命题 1.2（两个不交的核）加 FK 恰好就是它的标准证明。空隙情形（引理 2.4）等价于 λ₂(Ω₊) ≥ 4κ > Λ。

## 3. Pólya–Szegő 不等式（Schwarz 与 Steiner），含 H¹₀ 版本与条带版本

**内容。** 对非负函数做 Schwarz 或 Steiner 对称化，L² 范数不变，Dirichlet 积分不增。

**出处。**
- G. Pólya, G. Szegő, *Isoperimetric Inequalities in Mathematical Physics*, Ann. Math. Studies 27, Princeton 1951。
- Schwarz：Lieb–Loss Lemma 7.17；J. Brothers, W. Ziemer, J. Reine Angew. Math. 384 (1988) 153–179。
- Steiner：B. Kawohl, *Rearrangements and Convexity of Level Sets in PDE*, LNM 1150, Springer 1985, Ch. II（章节级）；F. Brock, A. Yu. Solynin, Trans. AMS 352 (2000) 1759–1796（极化法）；A. Baernstein II, *Symmetrization in Analysis*, CUP 2019（章节级）。

**用处与假设核对。**
- **引理 2.2(ii)。** 作用于 u ∈ C_c^∞(W)，W 开、有限测度。这是教科书情形。|u|* 是 Lipschitz 函数，支在紧集 K* ⊂ S_c 内，S_c 开（c 下半连续），所以 |u|* ∈ H¹₀(S_c)。✓
- **引理 4.5 第 3 步（x 方向 Steiner）。** 作用对象是 v ∈ H¹(T)，T = ℝ × (−τ, τ)，**不是 H¹₀**，墙上不要求为 0。
  - 教科书版本多是 H¹₀ 或 H¹(ℝⁿ)，所以这里需要一句局部化。
  - 用极化法即可：关于竖直线的极化只交换同一高度上的点，所以对任意水平带 ℝ × J，∫_{ℝ×J}|∇u^H|² = ∫_{ℝ×J}|∇u|²。
  - 再用迭代极化 L² 收敛到 Steiner 对称化（Brock–Solynin 2000），以及开集 ℝ × J 上能量的弱下半连续性，得到 ∫_T|∇v^♯|² ≤ ∫_T|∇v|²。
  - 每条水平线上 v(·, y) ∈ H¹(ℝ)，所以水平集测度有限。✓（建议在 proof 中写出这一句，见 解析稿 S2。）
- **引理 4.2、4.7（带轮廓的 PS）：** 见第 5 条。

## 4. 一个变量的单调重排（无边界条件）

**内容。** 在柱体 Ω′ × (0, τ) 上，对最后一个变量做单调递减重排，H¹ 范数保持可控，Dirichlet 积分不增；区间端点处不需要边界条件。

**出处。** Kawohl LNM 1150, Ch. II（章节级）；H. Berestycki, T. Lachand-Robert, Math. Nachr. 266 (2004) 3–19。

**用处与假设核对。** 用在引理 4.5 第 2 步：ℝ × (0, τ) 上的 y 向重排。
- 文献通常设 Ω′ 有界。这里 Ω′ = ℝ，但重排是逐纤维（固定 x）进行的，限制到 (−R, R) × (0, τ) 与重排可交换，令 R → ∞ 即可。✓
- 要求 v ≥ 0：proof 中应先把 v 换成 |v|（能量与迹的平方不变），见 解析稿 S2。
- 区间上的一维单调重排满足 PS，不需要边界条件（解析稿 §2 已说明：偶化之后每个非空非满的超水平集至少有一个内部边界点）。✓

## 5. 余面积公式与带等周轮廓的 Pólya–Szegő（Talenti 型能量不等式）

**内容。** 若对 a.e. t，Per({u > t}; G) ≥ I(|{u > t}|)，则 ∫_G|∇u|² ≥ ∫ I(s)²(u*′(s))² ds。证明只用余面积公式加 Cauchy–Schwarz。

**出处。**
- H. Federer, *Geometric Measure Theory*, Springer 1969, 3.2.12；L. C. Evans, R. F. Gariepy, *Measure Theory and Fine Properties of Functions*, CRC 1992（2015 修订版），§3.4（Lipschitz）与 §5.5（BV）。
- G. Talenti, Ann. Scuola Norm. Sup. Pisa (4) 3 (1976) 697–718；Brothers–Ziemer 1988（u* 的绝对连续性与平台处理）；S. Kesavan, *Symmetrization and Applications*, World Scientific 2006（章节级）。

**用处与假设核对。**
- **引理 4.2（颈部）。** u ∈ H¹₀(S_c) 是非负基态，零延拓后在 H¹(ℝ²) 中；轮廓 I_neck(s) ≥ √(4πs) > 0。✓ 不需要 S_c 正则。
- **引理 4.7（窗口）。** v ∈ H¹(T) 双单调，周长取相对于 T 的周长，余面积公式在开集 T 上成立。
  - 还需要 v* 在 (0, A₀] 上局部绝对连续。
  - 双 Steiner 水平集的相对周长有正下界：有满列时 ≥ 4τ（引理 4.6 第一项 ≥ 2·2τ），无满列时 = Per(U) ≥ √(4πs)。所以轮廓 ≥ min(√(4πs), 4τ) > 0。
  - 有了正下界，Brothers–Ziemer/Talenti 的论证就给出 v* 的局部绝对连续。✓（proof 中只写了"先截断"，建议补这一句，见 解析稿 S1。）

**说明。** 题目要求的"Talenti 比较"（PDE 比较原理版本）在本链中**没有用到**。用的只是上面的能量型不等式。

## 6. 等周不等式；Dido 与相对等周

**等周不等式（有限周长集）。**
- E. De Giorgi, Atti Accad. Naz. Lincei Mem. Cl. Sci. Fis. Mat. Nat. Sez. I (8) 5 (1958) 33–44；F. Maggi, *Sets of Finite Perimeter and Geometric Variational Problems*, CUP 2012, Ch. 14（章节级）；Evans–Gariepy §5.6。
- 用处：引理 4.2 的 √(4πs)；引理 4.6 第二项；引理 4.6 第一项中对 V ∪ K 的使用。对任意有限周长集成立，无界也可以。✓

**BV 的迹与粘合、切片。**
- 内容：Per(V ∪ K) = Per(V; H) + Per(K; Hᶜ) + ∫_{∂H}|tr V − tr K|；对 a.e. t，Per(E ∩ {y > t}) = Per(E; {y > t}) + H¹(E ∩ {y = t})。
- 出处：Evans–Gariepy §5.3–5.4（迹与延拓），Maggi 2012 Ch. 16、18（章节级）。
- 用处：引理 4.6 第一项（弦上的竖直落差恰为迹差）；引理 4.2 的切割恒等式（y_t 取在正测度的好切片集中，a.e. 成立）。✓

**Dido（定弦）。** 引理 4.6 第一项用"V ∪ 对侧优弓形 K"直接归结为等周不等式，对一切面积 A > 0 成立，不需要另引 Dido 定理。✓

**相对等周（半平面，√(2πs)）。** 本链中没有用到。它只出现在定理 A 的迹引理里，而定理 A 不在本链上。

## 7. Courant 节点定理，或替代

**本链不需要 Courant 定理。** 命题 1.2 只用到：
- (i) 本征函数在开集内连续（Weyl 引理，见第 10 条）；
- (ii) Hedberg 型引理：u ∈ H¹₀(U) ∩ C(U) ⟹ u⁺ ∈ H¹₀({u > 0})。proof 中已自证：取 min(u_n⁺, (u − ε)⁺)，利用 H¹₀ 的弱闭性，再令 ε → 0。所用的链式法则是 D. Gilbarg, N. S. Trudinger, *Elliptic PDE of Second Order*, Lemma 7.6；
- (iii) 正部、负部作试探函数，加上本征函数的正交性。

参考：R. Courant, D. Hilbert, *Methods of Mathematical Physics* I, Ch. VI §6（节点定理，仅作背景）。

**假设核对。** 偶子空间上 λ₁^{ev} 单重时 η₁ 不变号；λ₁^{ev} = λ₂^{ev} 时，在本征空间中任取两个正交本征函数。两种情况都得到两个不交的开核。✓（解析稿 S3 建议把这两种情况写开。）

## 8. Neumann 括号与 min-max

**出处。** M. Reed, B. Simon, *Methods of Modern Mathematical Physics* IV, Academic Press 1978, Thm XIII.1–XIII.2（min-max，形式版本）与 §XIII.15（Dirichlet–Neumann 括号）；Courant–Hilbert I, Ch. VI。

**用处与假设核对。** 命题 1.1：H¹₀(Ω) 经限制等距嵌入 V₊ ⊕ V₋（Ω ∩ L 是零测集），再用 min-max 比较。
- 只需要二次型是闭的、下有界。
- V± 的谱离散，因为 V± 经偶反射嵌入 H¹₀(S±)，而 S± 测度有限（第 9 条）。
- 不需要 Ω 正则或有界。✓

## 9. 有限测度开集上 H¹₀ ↪ L² 是紧嵌入

**内容与证明。** 设 |U| < ∞，U 可以无界。H¹₀(U) 零延拓后在 H¹(ℝ²) 中。
- 尾部：∫_{U∖B_R}u² ≤ ‖u‖²_{L⁴}|U ∖ B_R|^{1/2} ≤ C‖u‖²_{H¹}|U ∖ B_R|^{1/2} → 0，对单位球一致。
- B_R 内用 Rellich。

**出处。** Lieb–Loss Thm 8.5（二维 H¹ ⊂ L^q）与 Thm 8.6（Rellich–Kondrashov）。

**用处。** Ω、S±、W_i、S_c 的谱离散，基态存在（引理 4.2、4.4 需要）。所有这些集合测度都有限：|S_c| = ∫c ≤ 2a，|S±| = 2|Ω±|。✓

## 10. Weyl 引理 / 内部正则性

**出处。** H. Weyl, Duke Math. J. 7 (1940) 411–444；Gilbarg–Trudinger Thm 8.10（内部 H^k 正则性）加 Sobolev 嵌入。

**用处。** 命题 1.2：偶本征函数在 S₊ 内光滑，因而连续。只用内部正则性，与边界无关。✓

## 11. Sturm 理论、Prüfer 变换、Picone 恒等式

**出处。** H. Prüfer, Math. Ann. 95 (1926) 499–518；E. A. Coddington, N. Levinson, *Theory of ODE*, McGraw-Hill 1955, Ch. 8（章节级）；A. Zettl, *Sturm–Liouville Theory*, AMS Math. Surveys Monogr. 121 (2005)；M. Picone, Ann. Scuola Norm. Sup. Pisa 11 (1910) 1–144。

**用处与假设核对。**
- **判据 4.10，即证书 (S)。** 有限区间 [0, τ]，势 π²/c² ∈ [0, π²/w²] 有界（c 可以取 +∞，此时势为 0）。
  - 等价式 μ_h ≤ Λ₀ ⟺ θ(τ) ≥ θ_R 由基态代换自证；用到下确界可达（有界势、有限区间）。✓
- **颈部证书 N、证书 (T)。** 用的是 Riccati/Picone 超解判据，proof 中自证（我逐项展开核对过）。允许向上跳，作为测度处理。✓
- **定理 3.1，即证书 (B)。** 半直线问题。先把 c 截断到 [0, X]（μ 关于 X 的收敛来自紧支试探函数的稠密性），再加 ε′。
  - 之后 μ < Λ₀ ⟺ 0 端 Neumann 解在有限时间内到达 θ = π。这是标准振荡定理：有零点时 F·1_{[0,z]} 给出商 Λ₀，且它不是本征函数，所以下确界严格更小；无零点时用基态代换。
  - 两位复核者已核。✓

## 12. 其他基本工具

- Hölder、Cauchy–Schwarz：引理 4.3。
- Fatou：引理 2.2(i) 的下半连续性。
- Tonelli/Fubini：剖面、条带质量。
- 介值定理（S¹ 上的奇函数有零点）：定理 2.5。
- Sobolev 函数沿直线绝对连续（ACL）：W. P. Ziemer, *Weakly Differentiable Functions*, GTM 120 (1989), Thm 2.1.4；Evans–Gariepy §4.9。用于切片 Poincaré、帽内切片与横切片。
- H¹((−τ, τ); L²) ⊂ C([−τ, τ]; L²)：引理 4.5 中迹与重排可交换，引理 4.9 中 f(τ) 等于迹。

## 总表：每项假设是否满足

| 经典结果 | 要求 | 本链中的使用 | 是否满足 |
|---|---|---|---|
| FK（开集版） | 开、有限测度 | W_i、W_i⁺、半平面分支 | ✓ |
| FK（H¹ 版） | v ∈ H¹(ℝ²)，在有限测度集外为 0 | 帽反射、w* = 0 | ✓ |
| Steiner PS（C_c^∞ / H¹₀） | 开集、紧支 | 引理 2.2(ii) | ✓ |
| Steiner PS（条带 H¹） | 水平线上水平集测度有限 | 引理 4.5 第 3 步 | ✓（需写明局部化，S2） |
| 单调重排 | 非负、柱体、逐纤维 | 引理 4.5 第 2 步 | ✓（需先取 \|v\|，S2） |
| 带轮廓的 PS | 轮廓 > 0，u* 局部 AC | 引理 4.2、4.7 | ✓（引理 4.7 需补 AC 理由，S1） |
| 等周、BV 粘合与切片 | 有限周长 | 引理 4.2、4.6 | ✓ |
| min-max / 括号 | 闭的下有界二次型 | 命题 1.1 | ✓ |
| 紧嵌入 | 有限测度 | 全部谱的离散性 | ✓ |
| Weyl | 内部 | 命题 1.2 | ✓ |
| Prüfer / Sturm / Picone | 有界势或正则化 | (S)、(T)、N、(B) | ✓ |
| Courant | — | 不需要 | — |
| Talenti PDE 比较、相对等周 √(2πs)、Krahn–Szegő | — | 不直接使用 | — |
