"""Schematic of the spectral DCM generative model and the forward map F: A -> Sigma.

Four default-network nodes on an axial slice drive four coupled neuronal time series; the observed signal
is a (static) linear read-out (observation noise G_e = 0); the zero-lag covariance Sigma follows from the forward map.
Simplified redrawing after Fig. 2 of Friston et al. (2014, NeuroImage 99:533-547; CC BY 3.0), with a general
(asymmetric) effective-connectivity matrix A.  Time series, spectra and matrices are computed from the model
equations with illustrative parameter values (not results).
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, Rectangle
from scipy.linalg import solve_continuous_lyapunov

HERE = Path(__file__).resolve().parent
TITLE_Y = 4.70 + 0.06                   # y of the panel titles (the top of the panels is at 4.45 + DY)
DY = 0.0                                # title removed
W, H = 7.0, 5.10                       # canvas reduced to content + top margin (was 5.66); no top crop needed
INK, EDGE, FACE, GREY = "#222222", "#555555", "#f4f4f4", "#7f7f7f"
REGIONS = ["mPFC", "PCC", "lIPC", "rIPC"]
LABEL = {"lIPC": "L-IPC", "rIPC": "R-IPC"}      # display names (lIPC/rIPC = left/right inferior parietal cortex)
COLS = ["#E69F00", "#56B4E9", "#CC79A7", "#117733"]      # validated categorical set (light mode)
BETA, ALPHA = 0.5, 1.0
RAD = 0.11                              # curvature of the arrows in each reciprocal pair

plt.rcParams.update({"font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 9,
                     "axes.linewidth": 0.7, "pdf.fonttype": 42, "svg.fonttype": "none"})

# ---- illustrative model: A[i, j] = influence of region j on region i; mPFC drives PCC, lIPC and rIPC,
#      lIPC drives PCC, PCC drives rIPC, rIPC drives lIPC and feeds back to mPFC; no other couplings
#      (previous values: A[0,3]=0.25, A[1,0]=0.15, A[1,2]=0.60, A[2,0]=0.40, A[2,3]=0.35, A[3,0]=0.45, A[3,1]=0.60)
#      (set 2, kept as fallback: A[0,3]=0.30, A[1,0]=0.20, A[1,2]=0.50, A[2,0]=0.35, A[2,3]=0.40, A[3,0]=0.50, A[3,1]=0.45)
#      (set 3, all-positive off-diagonals: A[0,3]=0.20, A[1,0]=0.30, A[1,2]=0.70, A[2,0]=0.25, A[2,3]=0.50, A[3,0]=0.35, A[3,1]=0.40)
A = np.array([[-1.0, 0.00, 0.00, -0.30],
              [ 0.35, -1.0, 0.60, 0.00],
              [-0.25, 0.00, -1.0, 0.45],
              [ 0.40, -0.35, 0.00, -1.0]])
assert np.all(np.linalg.eigvals(A).real < 0)
Hm = np.eye(4)                           # H = I in this illustration
Sigma_w = np.array([[ 1.00, 0.45, 0.10, -0.10],
                    [ 0.45, 1.00, 0.30,  0.05],
                    [ 0.10, 0.30, 1.00,  0.35],
                    [-0.10, 0.05, 0.35,  1.00]])
assert np.all(np.linalg.eigvalsh(Sigma_w) > 0)

# ---- simulate  x' = A x + w  (power-law w),  y = H x  (G_e = 0)
rng = np.random.default_rng(7)
dt, T = 0.1, 400.0
n = int(T / dt)
f = np.fft.rfftfreq(n, dt); f[0] = f[1]
shape = f ** (-BETA / 2.0)
w_ind = np.fft.irfft(np.fft.rfft(rng.standard_normal((4, n)), axis=1) * shape, n=n, axis=1)
w_ind /= w_ind.std(axis=1, keepdims=True)
w = np.linalg.cholesky(Sigma_w) @ w_ind
x = np.zeros((4, n))
for k in range(n - 1):
    x[:, k + 1] = x[:, k] + dt * (A @ x[:, k] + w[:, k])
y = Hm @ x                                  # G_e = 0: no observation noise
t0, t1 = 2000, 2600                     # 60 s window after burn-in
sl = slice(t0, t1)

# ---- spectra and covariance from the same A (analytic)
def sigma_integral(A, beta, alpha, Sigma_w, Hm):
    """Sigma = H [alpha/(2 pi) int |w|^-beta R Sigma_w R^dagger dw] H^T."""
    from scipy.integrate import quad
    n_ = A.shape[0]; I_ = np.eye(n_)
    def entry(i, j, lo, hi, singular):
        def f(o):
            R = np.linalg.inv(1j * o * I_ - A)
            return (R @ Sigma_w @ R.conj().T)[i, j].real
        if singular:                                   # integrable |w|^-beta endpoint singularity at 0
            return quad(lambda o: f(o), lo, hi, weight="alg", wvar=(-beta, 0), limit=400)[0]
        return quad(lambda o: o ** (-beta) * f(o), lo, hi, limit=400)[0]
    B = np.zeros((n_, n_))
    for i in range(n_):
        for j in range(n_):
            B[i, j] = entry(i, j, 0.0, 1.0, True) + entry(i, j, 1.0, 400.0, False)
    return Hm @ (alpha / np.pi * B) @ Hm.T              # two-sided integral = 2 x one-sided, then 1/(2 pi)

Sigma = sigma_integral(A, BETA, ALPHA, Sigma_w, Hm)
assert np.allclose(Sigma, Sigma.T, atol=1e-8) and np.all(np.linalg.eigvalsh(Sigma) > 0)
# self-check: at beta = 0, A S + S A^T + alpha Sigma_w = 0
S0 = sigma_integral(A, 0.0, ALPHA, Sigma_w, np.eye(4))
assert np.allclose(S0, solve_continuous_lyapunov(A, -ALPHA * Sigma_w), atol=2e-3), "integral != Lyapunov at beta=0"

fig = plt.figure(figsize=(W, H))
bg = fig.add_axes([0, 0, 1, 1]); bg.set_xlim(0, W); bg.set_ylim(0, H); bg.axis("off")
def add_ax(x, y, w_, h_): return fig.add_axes([x / W, y / H, w_ / W, h_ / H])
def title(x, y, s): bg.text(x, y, s, ha="center", va="center", fontsize=11, color=INK)
def arrow(p, q, **kw): bg.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=10, lw=1.1, color=INK, **kw))

# ============================== top row: time domain
ty0, ty1 = 2.95 + DY, 4.45 + DY                   # vertical extent of the trace panels
# brain
meta = json.loads((HERE / "brain_axial.json").read_text())
img = plt.imread(HERE / "brain_axial.png")
bw_, bh_ = 1.62, 1.62 * img.shape[0] / img.shape[1]
bx0, by0 = 0.32, ty1 - bh_ + 0.12
bg.imshow(img, extent=[bx0, bx0 + bw_, by0, by0 + bh_], aspect="auto", zorder=1)
bg.set_xlim(0, W); bg.set_ylim(0, H)
def to_xy(mx, my):
    return (bx0 + (mx - meta["x_min"]) / (meta["x_max"] - meta["x_min"]) * bw_,
            by0 + (my - meta["y_min"]) / (meta["y_max"] - meta["y_min"]) * bh_)
P = [to_xy(meta["rois"][r]["x"], meta["rois"][r]["y"]) for r in REGIONS]
for i in range(4):                       # directed arrows j -> i, width proportional to |A[i, j]| (weak ones omitted; sign not shown)
    for j in range(4):
        if i != j and abs(A[i, j]) >= 0.15:
            rad = -2.0 * RAD if (i, j) == (2, 3) else RAD   # rIPC -> lIPC bows below the node line; the rest bow slightly
            bg.add_patch(FancyArrowPatch(P[j], P[i], arrowstyle="-|>", mutation_scale=8,
                                         lw=0.6 + 3.2 * abs(A[i, j]), color="#444444", zorder=3,
                                         connectionstyle=f"arc3,rad={rad}", shrinkA=6, shrinkB=6, alpha=0.9))
for (px, py), c, r in zip(P, COLS, REGIONS):
    bg.add_patch(Circle((px, py), 0.075, fc=c, ec="white", lw=1.0, zorder=5))
off = {"mPFC": (0, 0.17), "PCC": (-0.25, 0.03), "lIPC": (0.0, -0.17), "rIPC": (0.0, -0.17)}
for (px, py), r in zip(P, REGIONS):
    dx, dy = off[r]
    bg.text(px + dx, py + dy, LABEL.get(r, r), ha="center", va="center", fontsize=10, color=INK, zorder=6,
            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.75))

# the noise-spectrum equation, under the brain
cx = bx0 + bw_ / 2
blocks = [("neuronal noise spectra", r"$\mathbf{S}_w(\omega)=\alpha\,|\omega|^{-\beta}\,\boldsymbol{\Sigma}_w$")]
for q, (lab, eq_) in enumerate(blocks):
    y_ = (1.19 + DY + 1.15 / 2) + 0.41 - q * 0.44      # block (label to last equation, 0.82 in tall) centred on the matrices
    bg.text(cx, y_, lab, ha="center", va="center", fontsize=10, color=INK)
    bg.text(cx, y_ - 0.20, eq_, ha="center", va="center", fontsize=11, color=INK)
    bg.text(cx, y_ - 0.62, "observation noise spectra", ha="center", va="center", fontsize=10, color=INK)
    bg.text(cx, y_ - 0.82, r"$\mathbf{G}_e(\omega)=\mathbf{0}$", ha="center", va="center", fontsize=11, color=INK)

# time-series panels
px = [2.35, 3.95, 5.55]
pw = 1.15
data = [w[:, sl], x[:, sl], y[:, sl]]
names = ["Noise fluctuations", "Neuronal activity", "Observed BOLD"]
syms = [r"$\mathbf{w}(t)$", r"$\mathbf{x}(t)$", r"$\mathbf{y}(t)$"]
eqs = [r"$\mathbf{S}_w(\omega)=\alpha\,|\omega|^{-\beta}\,\boldsymbol{\Sigma}_w$", r"$\dot{\mathbf{x}}=\mathbf{A}\,\mathbf{x}+\mathbf{w}$",
       r"$\mathbf{y}=\mathbf{H}\,\mathbf{x}+\mathbf{e}$"]
for k in range(3):
    ax = add_ax(px[k], ty0, pw, ty1 - ty0)
    d = data[k] / data[k].std()
    lane = 2.4
    for r in range(4):
        ax.plot(np.arange(d.shape[1]) * dt, d[r] * 0.55 + (3 - r) * lane, color=COLS[r], lw=0.9)
    ax.set_xlim(0, d.shape[1] * dt); ax.set_ylim(-2.6, 3 * lane + 2.6)
    ax.axis("off")
    ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes, fc="none", ec=EDGE, lw=0.7, clip_on=False))
    title(px[k] + pw / 2, TITLE_Y + 0.13, names[k])
    lab_, fs_ = (syms[k], 11) if k == 0 else (eqs[k], 11)   # w(t) keeps its symbol; x and y show their model equations
    bg.text(px[k] + pw / 2, TITLE_Y - 0.10, lab_, ha="center", va="center", fontsize=fs_, color=INK)
for k in range(2):
    arrow((px[k] + pw + 0.07, (ty0 + ty1) / 2), (px[k + 1] - 0.07, (ty0 + ty1) / 2))

# title removed; now lives in the LaTeX caption (bolded)

# ============================== second row: one matrix under each trace panel
def chips(x0, y0, s_):
    """region-coloured chips along the left and top edge of a 4x4 square with lower-left corner (x0, y0)."""
    c = s_ / 4
    for r in range(4):
        bg.add_patch(Rectangle((x0 - 0.10, y0 + s_ - (r + 1) * c + 0.035), 0.06, c - 0.07, fc=COLS[r], ec="none"))
        bg.add_patch(Rectangle((x0 + r * c + 0.035, y0 + s_ + 0.035), c - 0.07, 0.06, fc=COLS[r], ec="none"))

ms = 1.15                                # matrix side (in)
yb = 1.19 + DY
mats = [Sigma_w, A, Sigma]
vmax = [1.0, 1.0, np.abs(Sigma).max()]
cmaps = ["RdBu_r", "RdBu_r", "RdBu_r"]
vmins = [-1.0, -1.0, -vmax[2]]
mnames = [r"$\boldsymbol{\Sigma}_w$ (correlated)", r"$\mathbf{A}$ (asymmetric)",
          r"$\boldsymbol{\Sigma}$ (symmetric)"]
for k in range(3):
    x0 = px[k] + pw / 2 - ms / 2
    ax = add_ax(x0, yb, ms, ms)
    ax.imshow(mats[k], cmap=cmaps[k], vmin=vmins[k], vmax=vmax[k], interpolation="nearest")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_linewidth(0.7)
    chips(x0, yb, ms)
    bg.text(px[k] + pw / 2, yb - 0.20, mnames[k], ha="center", va="center", fontsize=11, color=INK)
    arrow((px[k] + pw / 2, ty0 - 0.10), (px[k] + pw / 2, ty0 - 0.10 - 0.31))   # trace panel -> its matrix (same length as the horizontal arrows)

# ============================== third row: the forward map
bg.text(W / 2, 0.60,
        r"$F:\ \mathbf{A}\ \longmapsto\ \boldsymbol{\Sigma}(\mathbf{A},\beta)$",
        ha="center", va="center", fontsize=13, color=INK)

from matplotlib.transforms import Bbox
CROP_BOTTOM = 0.24                      # inches trimmed from the bottom of the canvas (nothing is drawn there)
CROP_TOP = 0.0                          # no top crop; canvas already sized to content
crop = Bbox.from_bounds(0.0, CROP_BOTTOM, W, H - CROP_BOTTOM - CROP_TOP)
fig.savefig(HERE / "spectral_dcm_schematic.pdf", dpi=300, bbox_inches=crop)
fig.savefig(HERE / "spectral_dcm_schematic.png", dpi=300, bbox_inches=crop)
# fig.savefig(HERE / "spectral_dcm_schematic.svg", dpi=300, bbox_inches=crop)
print("written")
