# Nomenclature

**Project:** A Lyapunov factorization of the map from causal to correlational brain connectivity at rest: intrinsic dynamics, temporal coloring, and spatial covariance  
**Last updated:** 2026-11-18

Single source of truth for mathematical notation, symbol definitions, and
abbreviations. Update this registry whenever
a symbol or term is introduced, removed, or given a different meaning.

**Scope tags** in the Notes columns: *Main* = main manuscript; *App A* =
dyadic adjugate derivation; *App B* = general-network covariance and
Lyapunov derivation; *App C* = systematic review. Some older row-level
scope tags have not yet been rechecked against this appendix mapping.

---

## 1. Typographic conventions

| Category | Rule | Examples |
|----------|------|----------|
| Matrices | Bold uppercase; bold Greek for Greek matrices | $\mathbf{A}$, $\mathbf{H}$, $\boldsymbol{\Sigma}$, $\boldsymbol{\Lambda}$ |
| Vectors and vector-valued processes | Bold lowercase | $\boldsymbol{x}(t)$, $\boldsymbol{y}(t)$, $\boldsymbol{w}(t)$, $\boldsymbol{e}(t)$, $\mathbf{e}_j$ |
| Scalars | Plain italic | $\alpha$, $\beta$, $\omega$, $a$, $f$, $g$ |
| Scalar integrals | Calligraphic | $\mathcal{J}$, $\mathcal{L}$ in App B |
| Sets and number systems | Calligraphic for matrix sets; blackboard bold for number systems | $\mathcal{U}$, $\mathcal{K}$; $\mathbb{R}$, $\mathbb{C}$ |
| Operators | Upright roman | $\operatorname{tr}$, $\operatorname{adj}$, $\operatorname{vec}$, $\operatorname{vech}$, $\operatorname{Re}$ |
| Transpose | Superscript $\top$ | $\mathbf{A}^{\top}$ |
| Conjugate transpose | Superscript $\dagger$ | $\mathbf{P}^{\dagger}$, $\mathbf{R}^{\dagger}$ |
| Complex conjugate | Overline | $\bar\lambda_l$ |
| Matrix entries | Matching lowercase indices | $(M_w)_{kl}$, $M_{kl}$, $A_{ij}$ |
| Spelling | American English | modeling, diagonalizable, colored noise |

---

## 2. Operators and decorations

| Symbol | Meaning | Notes |
|--------|---------|-------|
| $\operatorname{tr}(\mathbf{X})$ | Trace of $\mathbf{X}$ | Macro `\tr`; Main, App A, App B |
| $\det(\mathbf{X})$ | Determinant of $\mathbf{X}$ | Main dyad and App B |
| $\operatorname{adj}(\mathbf{X})$ | Adjugate of $\mathbf{X}$ | Macro `\adj`; for a $2\times2$ matrix, $\operatorname{adj}(\mathbf{X})=\operatorname{tr}(\mathbf{X})\mathbf{I}-\mathbf{X}$; Main, App B |
| $\operatorname{vec}(\mathbf{X})$ | Column-stacking vectorization | Macro `\vc`; App A, App B |
| $\operatorname{vech}(\mathbf{X})$ | Half-vectorization of a symmetric matrix | Macro `\vech`; collects the lower triangle and diagonal once; App A |
| $\operatorname{Re}z$ | Real part of a complex number | Written `\operatorname{Re}` or macro `\re`; App A, App B |
| $\operatorname{diag}(\cdot)$ | Diagonal matrix formed from its arguments | Main, App A |
| $\mathbb{E}[\cdot]$ | Expectation | Macro `\E`; App A |
| ${}^{\top}$ | Transpose | Used for real matrices and vectors throughout |
| ${}^{\dagger}$ | Conjugate transpose | Not a pseudoinverse; Main, App A, App B |
| $\bar z$ | Complex conjugate of $z$ | Most often $\bar\lambda_l$ |
| $\circ$ | Hadamard, or entrywise, matrix product | App A only |
| $\|\cdot\|$ | Spectral norm for matrices | App A convergence argument |
| $[x_{kl}]_{k,l=1}^{n}$ | Matrix whose $(k,l)$ entry is $x_{kl}$ | Main, App A |
| $(-\mathbf{A})^{-\beta}$ | Principal matrix power | Defined from the principal scalar powers of the eigenvalues of $-\mathbf{A}$; Main, App A |
| $\operatorname{Log}z$ | Principal complex logarithm | Used implicitly in principal powers; App B |
| $\Gamma(\cdot)$ | Euler gamma function | App B elementary integral |
| $B(p,q)$ | Euler beta function | $B(p,q)=\Gamma(p)\Gamma(q)/\Gamma(p+q)$; App B |

---

## 3. Symbols

### 3.1 Scalars and indices

| Symbol | Meaning | Notes |
|--------|---------|-------|
| $n$ | Number of brain regions or network nodes | $n\ge2$ in the formal setting; Main, App A |
| $t$ | Continuous time | Main, App A |
| $i$ | Imaginary unit, $i^2=-1$ | All mathematical sections |
| $j,k,l$ | Region, eigenmode, or matrix-entry indices | Meaning follows the local sum or entry |
| $s_j$ | Neuronal-noise amplitude for region $j$ | $s_j\ge0$; App B |
| $a$ | Half trace of the dyadic connectivity matrix | $a:=\operatorname{tr}(\mathbf{A})/2<0$ under stability; Main, App B |
| $f$ | First dyadic divided difference | $f=[(-\lambda_2)^{1-\beta}-(-\lambda_1)^{1-\beta}]/(\lambda_2-\lambda_1)$; real and negative; Main, App B |
| $g$ | Second dyadic divided difference | $g=[(-\lambda_1)^{-1-\beta}-(-\lambda_2)^{-1-\beta}]/(\lambda_2-\lambda_1)$; real and negative; Main, App B |
| $q_{kl}$ | Scalar integral for an eigenvalue pair | $q_{kl}=\int_{\mathbb R}|\omega|^{-\beta}/[(i\omega-\lambda_k)(-i\omega-\bar\lambda_l)]\,d\omega$; App A |
| $p,q$ | Arguments of the Euler beta function | Local to App B |
| $c$ | Generic complex parameter in the elementary integral | $\operatorname{Re}c<0$; App A |
| $K,C,\delta,m$ | Positive bounding constants | Local to convergence and holomorphy arguments; App A |
| $N$ | Frequency cutoff | Symmetric integration over $[-N,N]$; App A |

### 3.2 Vectors and vector-valued processes

| Symbol | Meaning | Notes |
|--------|---------|-------|
| $\boldsymbol{x}(t)$ | Neuronal state vector | $\in\mathbb{R}^n$; Main, App A |
| $\boldsymbol{y}(t)$ | Observed BOLD-signal vector | $\in\mathbb{R}^n$; Main, App A |
| $\boldsymbol{w}(t)$ | Neuronal-noise process | Stationary generalized input; Main, App A |
| $\boldsymbol{e}(t)$ | Observation-noise process | $\in\mathbb{R}^n$; set to zero spectrally through $\mathbf{G}_e\equiv\mathbf{0}$; Main, App A |
| $\mathbf{e}_j$ | $j$th standard basis vector | Used to isolate region-specific noise components; Main, App A |

### 3.3 Matrices and matrix-valued functions

| Symbol | Meaning | Notes |
|--------|---------|-------|
| $\mathbf{A}$ | Effective-connectivity matrix | Real $n\times n$ matrix; assumed Hurwitz stable; Main, App A, App B |
| $\mathbf{H}$ | Observation or hemodynamic-response matrix | Fixed, real, time-invariant, and frequency-independent; set to $\mathbf I$ in the main derivation so neuronal and observed covariances coincide |
| $\mathbf{I}$ | Identity matrix | Size inferred from context |
| $\boldsymbol{\Sigma}(\mathbf{A},\alpha,\beta,\boldsymbol\Sigma_w)$ | Stationary zero-lag neuronal covariance | Forward-map output; also the observed covariance under the main-text specialization $\mathbf H=\mathbf I$; abbreviated $\boldsymbol{\Sigma}(\mathbf{A},\beta)$ when $\alpha$ and $\boldsymbol\Sigma_w$ are fixed; Main, App A |
| $\mathbf D_{\boldsymbol\Sigma}$ | Diagonal matrix of marginal variances | $\operatorname{diag}(\Sigma_{11},\ldots,\Sigma_{nn})$; used to normalize covariance to Pearson correlation; Main |
| $\mathbf C$ | Pearson correlation matrix | $\mathbf D_{\boldsymbol\Sigma}^{-1/2}\boldsymbol\Sigma\mathbf D_{\boldsymbol\Sigma}^{-1/2}$ when every marginal variance is positive; Main |
| $\boldsymbol{\Sigma}_w$ | Spatial covariance factor of the multivariate Gaussian generalized neuronal noise | Real symmetric positive semidefinite; diagonal for independent regional inputs; identity for isotropic noise |
| $\boldsymbol{\Sigma}_x,\boldsymbol{\Sigma}_y$ | Neuronal and observed covariances | $\boldsymbol{\Sigma}_y=\mathbf H\boldsymbol{\Sigma}_x\mathbf H^\top$; Results and Appendices A/B abbreviate $\boldsymbol{\Sigma}_x$ as $\boldsymbol{\Sigma}$ |
| $\sigma_{jk}$ | Entry of the spatial neuronal-noise factor | $(\boldsymbol\Sigma_w)_{jk}=\sigma_{kj}$; App B |
| $\beta_{jk}$ | Exponent of a spatially correlated input pair | $(\beta_j+\beta_k)/2$, in $(-1,1)$; App B |
| $\mathbf Q_{jk}$ | Symmetric spatial-noise contribution of pair $(j,k)$ | $\sigma_{jj}\mathbf e_j\mathbf e_j^\top$ if $j=k$; $\sigma_{jk}(\mathbf e_j\mathbf e_k^\top+\mathbf e_k\mathbf e_j^\top)$ if $j<k$; off-diagonal terms may be indefinite; App B |
| $\mathbf E_{jk}$ | Matrix selecting entry $(j,k)$ | $\mathbf e_j\mathbf e_k^\top$; App B |
| $\mathbf C_{jk}$ | Pair-selection matrix in vectorized coordinates | Satisfies $\operatorname{vec}(\mathbf Q_{jk})=\mathbf C_{jk}\operatorname{vec}(\boldsymbol\Sigma_w)$ for symmetric $\boldsymbol\Sigma_w$; App B |
| $\mathbf D(\omega)$ | Regional spectral scaling matrix | $\operatorname{diag}(|\omega|^{-\beta_j/2})$; remains diagonal when $\boldsymbol\Sigma_w$ is spatially correlated; App B |
| $\mathbf K(\mathbf A)$ | Kronecker representation of the Lyapunov operator | $\mathbf I\otimes\mathbf A+\mathbf A\otimes\mathbf I$; App B |
| $\mathcal P_{\mathbf A}$ | Intrinsic dynamical-propagation factor | $-[\mathbf I\otimes\mathbf A+\mathbf A\otimes\mathbf I]^{-1}$; represents mean-reverting intrinsic dynamics; Main |
| $\mathcal T_{\beta,\mathbf A}$ | Connectivity-dependent temporal factor | $\eta_\beta[\mathbf I\otimes(-\mathbf A)^{-\beta}+(-\mathbf A)^{-\beta}\otimes\mathbf I]$; represents scale-free temporal coloring; Main |
| $\mathbf v_w$ | Vectorized spatial covariance factor | $\operatorname{vec}(\boldsymbol\Sigma_w)$; represents the spatial covariance structure of the multivariate Gaussian in main-text Eq. (3) |
| $\mathbf B$ | Adjugate of two-region connectivity | $\operatorname{adj}(\mathbf A)=\operatorname{tr}(\mathbf A)\mathbf I-\mathbf A$; App B |
| $\mathbf{S}_w(\omega)$ | Neuronal-noise cross-spectral density | Common-exponent form $\alpha|\omega|^{-\beta}\boldsymbol{\Sigma}_w$; regional form $\alpha\mathbf D(\omega)\boldsymbol\Sigma_w\mathbf D(\omega)$ in App B |
| $\mathbf{G}_e(\omega)$ | Observation-noise cross-spectral density | Set to $\mathbf{0}$ throughout; Main, App A |
| $\mathbf{G}_y(\omega)$ | Cross-spectral density of the observed signal | Main, App A |
| $\mathbf{R}(\omega;\mathbf{A})$ | Resolvent of $\mathbf{A}$ | $(i\omega\mathbf{I}-\mathbf{A})^{-1}$; the matrix argument is omitted when clear; App A, App B |
| $\mathbf{P}$ | Matrix of right eigenvectors of $\mathbf{A}$ | Complex and invertible for diagonalizable $\mathbf{A}$; Main, App A |
| $\boldsymbol{\Lambda}$ | Diagonal eigenvalue matrix | $\operatorname{diag}(\lambda_1,\ldots,\lambda_n)$; Main, App A |
| $\mathbf{P}^{-\dagger}$ | Conjugate transpose of $\mathbf{P}^{-1}$ | Defined as $(\mathbf{P}^{-1})^{\dagger}$; Main, App A |
| $\mathbf{M}_w$ | Noise covariance in eigenvector coordinates | $\mathbf{P}^{-1}\boldsymbol{\Sigma}_w\mathbf{P}^{-\dagger}$; Main, App A |
| $\mathbf{M}$ | Isotropic-noise version of $\mathbf{M}_w$ | $\mathbf{P}^{-1}\mathbf{P}^{-\dagger}$; Main, App A |
| $\mathbf{M}^{(jk)}$ | Contribution of pair $(j,k)$ in eigenvector coordinates | $\mathbf{P}^{-1}\mathbf Q_{jk}\mathbf{P}^{-\dagger}$; App B |
| $\tilde{\boldsymbol{\Sigma}}$ | Covariance in eigenvector coordinates | $\mathbf{P}^{-1}\boldsymbol{\Sigma}\mathbf{P}^{-\dagger}$; App A |
| $\boldsymbol{\Sigma}_{xw}$ | Spectral state–noise cross term | Symmetric principal value when $\beta\le0$; formal zero-lag cross-covariance; equals $\eta_\beta(-\mathbf A)^{-\beta}\boldsymbol{\Sigma}_w$ for every stable $\mathbf A$ |
| $\mathbf X$ | Generic matrix or difference between Lyapunov solutions | Generic $2\times2$ matrix in the adjugate identity (App B); solution difference in the uniqueness proof (App A) |
| $\mathbf N(\omega)$ | Numerator of $\mathbf R(\omega)\mathbf R(\omega)^\dagger$ | $\omega^2\mathbf I+\mathbf B\mathbf B^\top+i\omega(\mathbf B-\mathbf B^\top)$; App B |
| $\mathbf N_w(\omega)$ | Numerator of $\mathbf R(\omega)\boldsymbol\Sigma_w\mathbf R(\omega)^\dagger$ | $\omega^2\boldsymbol\Sigma_w+\mathbf B\boldsymbol\Sigma_w\mathbf B^\top+i\omega(\mathbf B\boldsymbol\Sigma_w-\boldsymbol\Sigma_w\mathbf B^\top)$, where $\mathbf B=\operatorname{adj}(\mathbf A)$; App A |
| $\mathbf Q_\beta(\mathbf A,\boldsymbol\Sigma_w)$ | Connectivity-dependent colored-noise forcing in the Lyapunov equation | $\eta_\beta[(-\mathbf A)^{-\beta}\boldsymbol\Sigma_w+\boldsymbol\Sigma_w[(-\mathbf A)^{-\beta}]^\top]$; App B |
| $\mathbf B_\beta$ | Power-law spreading matrix | $(-\mathbf A)^{-\beta}/\cos(\pi\beta/2)$; used in the Discussion to compare the spectral forcing with the time-domain colored-noise Lyapunov form of Gradziuk et al.; Main |

---

## 4. Greek symbols and scalar functions

| Symbol | Meaning | Notes |
|--------|---------|-------|
| $\alpha$ | Neuronal-noise intensity | $\alpha>0$; Main, App A, App B |
| $\beta$ | Common power-law exponent of neuronal noise | $-1<\beta<1$; $\beta=0$ is white noise; Main, App A, App B |
| $\beta_j$ | Region-specific power-law exponent | $-1<\beta_j<1$; App B |
| $\omega$ | Angular-frequency variable | Real; Main, App A, App B |
| $\lambda_k$ | $k$th eigenvalue of $\mathbf{A}$ | $\operatorname{Re}\lambda_k<0$ under Hurwitz stability; Main, App A, App B |
| $\lambda,\mu$ | Generic eigenvalue pair in the scalar-integral derivation | $\lambda=\lambda_k$, $\mu=\bar\lambda_l$; App A |
| $\eta_\beta$ | Colored-noise scale factor | $\eta_\beta:=\alpha/[2\cos(\pi\beta/2)]$; Main, App A |
| $\Delta(\omega)$ | Dyadic characteristic determinant | $\det(i\omega\mathbf{I}-\mathbf{A})=(i\omega-\lambda_1)(i\omega-\lambda_2)$; App B |
| $\mathcal{J}$ | Dyadic scalar integral multiplying $\mathbf{I}$ | $\int |\omega|^{-\beta}\omega^2/|\Delta(\omega)|^2\,d\omega$; App B |
| $\mathcal{L}$ | Dyadic scalar integral multiplying $\operatorname{adj}(\mathbf{A})\operatorname{adj}(\mathbf{A})^\top$ | $\int |\omega|^{-\beta}/|\Delta(\omega)|^2\,d\omega$; App B |
| $\zeta$ | Beta-function parameter | $\zeta=(1+\beta)/2$; local to App B |
| $\varphi(\lambda)$ | Repeated-eigenvalue auxiliary function | $(-\lambda)^{1-\beta}$; local to App B |
| $\psi(\lambda)$ | Repeated-eigenvalue auxiliary function | $(-\lambda)^{-1-\beta}$; local to App B |
| $\kappa$ | White-noise dyadic scale factor | $-\alpha/[2\operatorname{tr}(\mathbf{A})\det\mathbf{A}]$; local to App B |
| $\delta_{kl}$ | Kronecker delta | Used to select diagonal entries; App A |

---

## 5. Sets, spaces, and maps

| Symbol | Meaning | Notes |
|--------|---------|-------|
| $\mathbb{R}$ | Real numbers | Macro `\R` |
| $\mathbb{C}$ | Complex numbers | Used for eigenvectors, eigenvalues, and principal powers |
| $\mathcal{U}$ | Set of real Hurwitz-stable $n\times n$ matrices | Open domain of the forward map; App A |
| $F$ | Forward map from effective to functional connectivity | Written schematically as $F:\mathbf{A}\mapsto\boldsymbol{\Sigma}(\mathbf{A},\beta)$; vectorized form $F:\operatorname{vec}(\mathbf{A})\in\mathbb{R}^{n^2}\mapsto\operatorname{vech}\boldsymbol{\Sigma}\in\mathbb{R}^{n(n+1)/2}$; Main, App A |

---

## 6. Abbreviations

| Abbreviation | Expansion | Notes |
|--------------|-----------|-------|
| DCM | Dynamic causal modeling | Spectral DCM is the resting-state frequency-domain formulation used here |
| BOLD | Blood-oxygen-level-dependent | Observed fMRI signal modeled by $\boldsymbol{y}(t)$ |
| fMRI | Functional magnetic resonance imaging | Main and App C |
| EC | Effective connectivity | Directed/model-based connectivity; matrix $\mathbf{A}$ |
| FC | Functional connectivity | Correlational connectivity; covariance $\boldsymbol{\Sigma}$ |
| PRISMA | Preferred Reporting Items for Systematic Reviews and Meta-Analyses | Flow diagram in App C |
| mPFC | Medial prefrontal cortex | Figure 1 |
| PCC | Posterior cingulate cortex | Figure 1 |
| L-IPC / R-IPC | Left/right inferior parietal cortex | Figure 1 |
| IGD | Internet gaming disorder | Figure 2 |
| TBI | Traumatic brain injury | Figure 2 |
| OCD | Obsessive-compulsive disorder | Figure 2 |
