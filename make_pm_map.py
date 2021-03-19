import sys
sys.path.append("../../Source/pyflct")
import os
from pathlib import Path
import numpy as np
import pyflct
from astropy.io import fits
from astropy.wcs import WCS

os.environ["PATH"] = ":".join((os.environ["PATH"], "/Users/will/.local/bin"))

try:
    DATADIR, FN1, FN2, OUT_PREFIX = sys.argv[1:5]
except:
    sys.exit(f"Usage: {sys.argv[0]} DATADIR FN1 FN2 OUT_PREFIX [SIGMA]")

try:
    SIGMA = float(sys.argv[5])
except:
    # Optional input argument
    SIGMA = 10.0



data_path = Path(DATADIR)

h1, = fits.open(data_path / f"{FN1}.fits")
h2, = fits.open(data_path / f"{FN2}.fits")

# Remove NaNs in the input data
data1 = h1.data.T
data2 = h2.data.T
goodpixels = np.isfinite(data1*data2)
data1[~goodpixels] = 0.0
data2[~goodpixels] = 0.0
print("Number of goodpixels =", goodpixels.sum())

# Designed for images that have been high-pass filtered and normalized
data_scale = 2.0
sigma = SIGMA
vx, vy, vm = pyflct.flct(data1, data2,
                         deltat=1.0, deltas=1.0, sigma=sigma,
                         #thresh=0.1/data_scale,
)
#vx[vm==0] = np.nan
#vy[vm==0] = np.nan

fits.PrimaryHDU(
    header=h1.header,
    data=vx.T,
).writeto(data_path / f"{OUT_PREFIX}_vx_sig{int(sigma):02d}.fits", overwrite=True)

fits.PrimaryHDU(
    header=h1.header,
    data=vy.T,
).writeto(data_path / f"{OUT_PREFIX}_vy_sig{int(sigma):02d}.fits", overwrite=True)
