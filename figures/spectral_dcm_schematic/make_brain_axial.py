"""Brain asset for the spectral-DCM schematic: an axial slice of the MNI152 T1 template (nilearn, bundled data).

Anterior is up and the left hemisphere is on the LEFT (neurological convention).
Needs nilearn (run once; spectral_dcm_schematic.py itself only needs the PNG + JSON written here).
Writes brain_axial.png (RGBA, transparent outside the brain) and brain_axial.json
(extents in MNI mm, plus approximate default-network node centres in the same units; schematic positions).
"""
import json
import pathlib

import numpy as np
from nilearn import datasets
from PIL import Image
from scipy import ndimage

HERE = pathlib.Path(__file__).resolve().parent
Z_SLICE = 28                           # MNI z of the displayed slice (contains all four nodes)
ROIS = {"mPFC": (0, 38, 28), "PCC": (0, -40, 28),
        "lIPC": (-48, -66, 28), "rIPC": (50, -62, 28)}   # approximate, schematic default-network node positions (MNI mm), not fitted coordinates

t1 = datasets.load_mni152_template(resolution=1)
mask = datasets.load_mni152_brain_mask(resolution=1)
aff = t1.affine
iz = int(round(np.linalg.inv(aff).dot([0, 0, Z_SLICE, 1])[2]))
img = t1.get_fdata()[:, :, iz]           # [x, y]
msk = mask.get_fdata()[:, :, iz] > 0.5
msk = ndimage.binary_closing(msk, iterations=2)
msk = ndimage.binary_fill_holes(msk)

xs = np.nonzero(msk.any(1))[0]
ys = np.nonzero(msk.any(0))[0]
pad = 2
x0, x1 = xs.min() - pad, xs.max() + pad
y0, y1 = ys.min() - pad, ys.max() + pad
img, msk = img[x0:x1 + 1, y0:y1 + 1], msk[x0:x1 + 1, y0:y1 + 1]

UP = 6
hi = np.percentile(img[msk], 99.0)
g = np.clip(img / hi, 0, 1)
g = ndimage.zoom(g, UP, order=3)
m_up = ndimage.zoom(ndimage.gaussian_filter(msk.astype(float), 0.7), UP, order=3)
g = 0.30 + 0.68 * np.clip(g, 0, 1)       # lighten so coloured nodes and arrows read on top
alpha = np.clip((m_up - 0.35) / 0.30, 0, 1)
rgba = np.dstack([g, g, g, alpha])
rgba = np.transpose(rgba, (1, 0, 2))[::-1, :, :]   # rows = y (top = anterior), columns = x (left = -x)
Image.fromarray((rgba * 255).astype(np.uint8), "RGBA").save(HERE / "brain_axial.png")

x_mni = lambda ix: aff[0, 3] + ix * aff[0, 0]
y_mni = lambda iy: aff[1, 3] + iy * aff[1, 1]
meta = {"z_slice": Z_SLICE,
        "x_min": x_mni(x0 - 0.5), "x_max": x_mni(x1 + 0.5),
        "y_min": y_mni(y0 - 0.5), "y_max": y_mni(y1 + 0.5),
        "rois": {k: {"x": v[0], "y": v[1], "z": v[2]} for k, v in ROIS.items()}}
(HERE / "brain_axial.json").write_text(json.dumps(meta, indent=2))
print(json.dumps(meta, indent=2))
