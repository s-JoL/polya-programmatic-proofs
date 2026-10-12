# All-phase Q2/Q1 inverse-compression lemma

s-JoL · 2026-10-11 · Research preprint

The polynomial core is formalized in Lean. This exposition gives the full finite-matrix return using classical Fourier and Schur-complement arguments. No priority claim is made.

## 1. Deriving the wrapped symbol and every contraction

Let z=exp(it), c=(z+z^{-1})/2, s=1+z, and h=|s|²=2+2c. In coefficient order (vertex, midpoint), the map into one element is

    B = [[1,0],[0,1],[z,0]].

Substitution in the supplied raw element matrices gives, without any phase convention inferred from the parent,

    3 B* d B = D = [[14+2c, -8 conjugate(s)],[-8s,16]],
    15 B* m B = M = [[4-c, conjugate(s)],[s,8]].

Indeed, the vertex mass numerator before division by 30 is 8-z-z^{-1}, the vertex stiffness numerator before division by 3 is 14+z+z^{-1}, and the off-diagonal entries are respectively 2(1+z^{-1})/30 and -8(1+z^{-1})/3. Thus the supplied corrected blocks are exactly the assembled unit-cell symbols.

For j=(1,s/2) and g=(-conjugate(s)/2,1), j*g=0 and both squared norms are

    N = 1+h/4 = (3+c)/2 > 0.

The diagonal contractions, expanded before simplifying, are

    p = j*Mj = (4-c)+h+2h = 5(2+c),
    q = g*Mg = (4-c)h/4-h+8 = (16-c-c²)/2,
    r = j*Dj = 14+2c-8h+4h = 6(1-c),
    v = g*Dg = (14+2c)h/4+8h+16 = c²+24c+39.

The cross contractions are obtained directly from the two coordinates:

    a = j*Mg = conjugate(s)[-(4-c)/2+1+(8-h/2)/2]
      = (5/2) conjugate(s),
    b = j*Dg = conjugate(s)[-(14+2c)/2-8+(4h+16)/2]
      = -3(1-c) conjugate(s).

These identities hold also at s=0: no division by s occurs. The incorrect coefficient 8c in v is rejected. At c=1 the true v is 64, whereas c²+8c+39 would give 48; this is already an exact falsifier of that erroneous formula.

For all c in [-1,1], M11=4-c>=3 and det M=10(3-c)>=20, so M is positive definite. Also D11>=12 and det D=96(1-c)>=0, so D is positive semidefinite. These are Hermitian quadratic forms in the coefficient Euclidean inner product.

## 2. Tensor order, complete complement, and cross factors

Use y-major tensor coefficient order (vertex-y/vertex-x, vertex-y/midpoint-x, midpoint-y/vertex-x, midpoint-y/midpoint-x). Then the physical unit-cell pencil is

    A_E = d_y tensor m_x + m_y tensor d_x - E m_y tensor m_x,
    T_E := 225 A_E = 5(D_y tensor M_x+M_y tensor D_x)-E M_y tensor M_x.

An x-major convention is obtained by a unitary tensor permutation and changes none of the conclusions. The geometric Q1 prolongation is j_y tensor j_x, including the center value s_y s_x/4.

The four products of j and g are pairwise Euclidean orthogonal, all with squared norm N_y N_x>0. Consequently

    e1=j_y tensor g_x, e2=g_y tensor j_x, e3=g_y tensor g_x

are a complete basis of the three-dimensional orthogonal complement. These are not normalized. Positivity of their Gram restriction C_E=(e_i* T_E e_j) is equivalent to positivity on the actual complement.

With p,q,r,v,a,b as above and axis subscripts, every independent entry is

    C11 = 5(r_y q_x+p_y v_x)-E p_y q_x,
    C22 = 5(v_y p_x+q_y r_x)-E q_y p_x,
    C33 = 5(v_y q_x+q_y v_x)-E q_y q_x,
    C12 = 5(b_y conjugate(a_x)+a_y conjugate(b_x))-E a_y conjugate(a_x),
    C13 = 5(b_y q_x+a_y v_x)-E a_y q_x,
    C23 = 5(v_y a_x+q_y b_x)-E q_y a_x.

The remaining entries are their conjugates. In particular, putting x=c_x, y=c_y and E=8,

    C12 = conjugate(s_y) s_x F12, F12=(25/2)[3(x+y)-10],
    C13 = conjugate(s_y) F13, F13=(15y-35)q_x+(25/2)v_x,
    C23 = conjugate(s_x) F23, F23=(25/2)v_y+(15x-35)q_y.

The cross term in the determinant is real because

    C12 C23 conjugate(C13) = h_x h_y F12 F23 F13.

Writing A0=C11, B0=C22, C0=C33 at E=8, the leading minors are exactly

    Delta1 = A0,
    Delta2 = A0 B0 - h_x h_y F12²,
    Delta3 = A0 B0 C0 + 2 h_x h_y F12 F13 F23
             - A0 h_x F23² - B0 h_y F13² - C0 h_x h_y F12².

The independently implemented permutation determinant of the direct 4-by-4 tensor restriction equals this Hermitian determinant formula in the entire Laurent polynomial ring. Thus no modulus or cross-factor was dropped.

## 3. Exact all-phase positivity certificate

The checker derives D and M by matrix multiplication from the original three-node integer matrices, forms the four-dimensional tensor pencil, constructs the complement columns, and computes its leading principal minors. It uses rational Laurent polynomials in z,w; adjoint is exponent inversion. It verifies each 1D contraction and all six entries against a second contraction-level formula.

Each minor is invariant under z -> z^{-1} and w -> w^{-1} separately. Reduce each symmetric Laurent term using z^n+z^{-n}=2 T_n(x) and the corresponding y identity, with the recurrence 2T_{n+1}=2x(2T_n)-2T_{n-1}. Substitution of x=(z+z^{-1})/2 and y=(w+w^{-1})/2 back into every cosine polynomial reproduces the original Laurent determinant exactly. This specifically checks that the conversion to cosines has not hidden a phase dependence.

For Delta(2t-1,2u-1)=sum a_ij t^i u^j, its Bernstein coefficients at bidegree (n,m) are

    b_kl = sum_{i<=k,j<=l} a_ij binom(k,i)/binom(n,i)
                                  binom(l,j)/binom(m,j).

The identity t^i=sum_{k>=i}[binom(k,i)/binom(n,i)] B_k^n(t) follows directly from binom(n,k)binom(k,i)=binom(n,i)binom(n-i,k-i) and the binomial theorem. Tensoring this identity proves the conversion formula. Independently expanding sum b_kl B_k^n(t) B_l^m(u) reproduces the transformed polynomial coefficient by coefficient. The actual degrees and exact minimum coefficients are:

| Minor | Bidegree | Minimum Bernstein coefficient |
|---|---|---|
| Delta1 | (2,1) | 240 |
| Delta2 | (3,3) | 313600 |
| Delta3 | (5,5) | 240844800 |

All coefficients, including all 36 determinant coefficients, are retained in full in `certificates/torus.json`, together with every cosine power coefficient. For example the directly derived first minor is

    Delta1 = 1550+415y+1225x+635xy+75x²+60x²y,

with Bernstein coefficient rows (560,240), (1120,1830), (1740,3960). Bernstein basis functions are nonnegative and their sum is one on the closed square [0,1]². A strictly positive minimum therefore proves strict positivity everywhere, including all corners and edges. This is not a sampling argument or an extrapolation from numerical success.

By Sylvester's criterion the Hermitian C_8 is positive definite for every real pair of phases. As M_x,M_y are positive definite, their tensor product is positive definite (diagonalize the two positive matrices unitarily and multiply their positive eigenvalues). For E<=8,

    C_E = C_8+(8-E) V*(M_y tensor M_x)V,

where V has columns e1,e2,e3. Hence C_E>=C_8>0. This proves the requested interval 0<E<=8, without requiring full A_E itself to be positive or invertible.

## 4. Inertia and inverse compression: explicit congruence

This step is finite-dimensional and does not use a counting approximation. Let W=range J, with J full column rank, and assume A is Hermitian and positive on W-perp. Choose a Euclidean unitary U=[U0,U1] with U0 spanning W and write

    U* A U = [[H,K*],[K,C]], C>0,
    S = H-K* C^{-1}K.

For the invertible triangular matrix

    P = [[I,0],[-C^{-1}K,I]],

direct multiplication gives

    P* (U* A U) P = diag(S,C).

An invertible change of coordinates carries maximal negative subspaces and kernels bijectively, so this congruence preserves inertia and nullity. Thus n_-(A)=n_-(S), and A is invertible exactly when S is invertible. Solving the block linear equations for a right-hand side (v,0) shows that the top-left block of (U* A U)^{-1} is S^{-1}.

Write J=U0 R, where R is square and invertible. Therefore

    J* A^{-1} J = R* S^{-1} R.

Congruence by R preserves inertia. A nonzero real eigenvalue lambda of S has inverse 1/lambda with the same sign; hence

    n_-(J* A^{-1}J)=n_-(S^{-1})=n_-(S)=n_-(A),
    nullity(J* A^{-1}J)=0.

This argument also takes care of nonunit j/g norms. Identifying an inverse compression with the inverse of the Q1 restriction J*AJ would be wrong; the Schur correction is essential.

## 5. Finite twisted torus: assembly and full multiplicities

Let L>=1 be an integer, with L unit cells in each direction, no deleted cells, and boundary twists exp(i theta_x), exp(i theta_y). A species coefficient satisfies a_{n+L}=exp(i theta)a_n. The unitary one-dimensional Fourier vectors are

    f_a(n)=L^{-1/2} exp(i t_a n), t_a=(theta+2 pi a)/L,
    a=0,...,L-1.

Their orthonormality follows by summing the geometric series exp(2 pi i(a-b)n/L). There are exactly L vectors. Each cell has local vector (a_n,b_n,a_{n+1}); on a Fourier mode it is the common cell phase times B(t_a)(a,b). Summing local element quadratic forms over L cells cancels the normalization 1/L. Hence the assembled matrices have exactly the wrapped symbols above, with no extra L or multiplicity factor. This remains true for L=1 even though both endpoints identify through the twist.

Tensoring the two Fourier transforms and retaining both Q2 species in each direction yields exactly L² four-dimensional blocks, hence all 4L² fine coefficients. The assembled tensor FEM forms are the sums of the two derivative-mass terms and the mass product, so their blocks are precisely A_E(t_a,t_b). Applying the same Fourier transform to the L² coarse vertex values gives the column j_y(t_b) tensor j_x(t_a) in each block. Its four entries are (1,s_x/2,s_y/2,s_y s_x/4), exactly the bilinear interpolation values at the fine vertices, edges and center. The vertex component is one, so every block column is nonzero and the full J has rank L².

Consequently W-perp is the orthogonal direct sum of the three-dimensional complements already proved positive. Every negative eigenvalue of the full Hermitian pencil is accounted for by the direct sum of the L² blocks. Equal eigenvalues in distinct blocks retain their actual multiplicities; a conjugate phase is neither discarded nor spuriously duplicated. A single block has at most one nonpositive eigenvalue because its complement block C is positive. Applying Section 4 to the full assembled A and all coarse columns proves the claimed exact strict index equality and nullity zero whenever the full A is invertible.

For real periodic/antiperiodic data the matrix may also be regarded as a real symmetric matrix; complexification leaves its eigenvalue multiplicities unchanged. Fourier conjugate pairing is therefore a coordinate representation, not a doubling of the count.

## 6. Singular scope and discarded errors

The complement positivity does not need invertibility. At a singular Fourier fiber, the same congruence reduces that fiber to diag(s,C) with scalar s=0 and C>0. Its nullity is exactly one and its negative index is zero. An acoustic equality is consequently a zero mode; it does not introduce an extra strict negative eigenvalue in that fiber. Several distinct singular fibers contribute one kernel dimension each, solely because they are distinct orthogonal blocks. The full A is singular if and only if at least one block is singular. At such E the ordinary A^{-1} and its compression are undefined; no claim about a Moore-Penrose inverse, limiting compression, or closed-count equality has been established here.

The coefficient 24c in v is essential: the variant c²+8c+39 is wrong (Section 1), and the positivity certificate is computed from the exact full determinant, not from that variant.
