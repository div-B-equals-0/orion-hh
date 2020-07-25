import sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

try:
    DATADIR = sys.argv[1]
    SUFFIX = sys.argv[2]
    OUTDIR = sys.argv[3]
except IndexError:
    sys.exit(f"Usage: {sys.argv[0]} DATADIR SUFFIX OUTDIR")

infile = f"muse-hr-data-wavsec0{SUFFIX}-cont-sub.fits"
outfile = f"o_ii-v1-sum{SUFFIX}.fits"
hdu = fits.open(f"{DATADIR}/{infile}")["DATA"]
w = WCS(hdu)
nwav, ny, nx = hdu.data.shape


# Wavelength sections of O II V1 multiplet
v1_sections = [
    [49, 71], [79, 83], [90, 102]
]

im = np.zeros((ny, nx))
for i1, i2 in v1_sections:
    im += np.sum(hdu.data[i1:i2, :, :], axis=0)

fits.PrimaryHDU(header=w.celestial.to_header(), data=im).writeto(
        f"{OUTDIR}/{outfile}", overwrite=True)
print(f"Written {outfile}")
