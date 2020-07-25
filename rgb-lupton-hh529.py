STRETCH=0.9
Q=2.5
WINDOW=[[1354, 1369, 760, 600]]
TAB=[["final502-radec", 130, 750], ["final658-radec", 200, 2200], ["final673-radec", 70, 500]]
DIR="/Users/will/Work/BobPC/2002"
SUFFIX="1996"
HDU=0
from pathlib import Path
import matplotlib.pyplot as plt
from astropy.visualization import make_lupton_rgb
from astropy.io import fits
from astropy.wcs import WCS


datadir = Path(DIR)
# Unpack the channel info from the table
[rf, r1, r2], [gf, g1, g2], [bf, b1, b2] = TAB

def load_and_scale_image(fn, v1, v2, ihdu=HDU):
    data = fits.open(datadir / f"{fn}.fits")[ihdu].data
    return (data - v1) / (v2 - v1)

w = WCS(fits.open(datadir / f"{rf}.fits")[HDU].header)

image_r = load_and_scale_image(rf, r1, r2)
image_g = load_and_scale_image(gf, g1, g2)
image_b = load_and_scale_image(bf, b1, b2)
image = make_lupton_rgb(image_r, image_g, image_b, stretch=STRETCH, Q=Q)

x0, y0, dx, dy = WINDOW


figfile = f"rgb-lupton-hh529-{SUFFIX}.jpg"
fig, ax = plt.subplots(figsize=(6.2, 5), subplot_kw=dict(projection=w))
ax.imshow(image)
ax.set(
    xlabel="RA (J2000)",
    ylabel="Dec (J2000)",
    xlim=[x0-dx/2, x0+dx/2],
    ylim=[y0-dy/2, y0+dy/2],
)
fig.tight_layout(rect=[0.15, 0.07, 1.0, 1.0])
fig.savefig(figfile, dpi=200)
