import sys
from pathlib import Path
import numpy as np
from astropy.io import fits

datadir = Path("~/Dropbox/Orion-HH-data/MUSE").expanduser()

try:
    fileroot = sys.argv[1]
    threshold = float(sys.argv[2])
    STEEPNESS = float(sys.argv[3])
    SUFFIX = sys.argv[4]
except IndexError:
    sys.exit(f'Usage: {sys.argv[0]} FILEROOT THRESHOLD STEEPNESS SUFFIX')


refhdu = fits.open(datadir / "linesum-O_III-5007-bin001.fits")["SCALED"]
outim = np.empty_like(refhdu.data)
nlist = [1, 2, 4, 8, 16, 32, 64, 128, 256]
for n in reversed(nlist):
    hdu = fits.open(datadir / f"{fileroot}-bin{n:03d}.fits")["SCALED"]
    mask = refhdu.data >= threshold/n**STEEPNESS
    outim[mask] = hdu.data[mask]

fits.PrimaryHDU(
    header=hdu.header,
    data=outim,
).writeto(
    datadir / f"{fileroot}-multibin{SUFFIX}.fits",
    overwrite=True,
)
