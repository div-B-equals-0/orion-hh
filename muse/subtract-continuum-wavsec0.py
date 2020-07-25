import sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from numpy.polynomial import Chebyshev as T
import itertools

try:
    DATADIR = sys.argv[1]
    SUFFIX = sys.argv[2]
    OUTDIR = sys.argv[3]
except IndexError:
    sys.exit(f"Usage: {sys.argv[0]} DATADIR SUFFIX OUTDIR")

infile = f"muse-hr-data-wavsec0{SUFFIX}.fits"
hdu = fits.open(f"{DATADIR}/{infile}")["DATA"]
w = WCS(hdu)
nwav, ny, nx = hdu.data.shape
wavpix = np.arange(nwav)

# Two pairs of adjacent sections for the true continuum

# Wavelength sections of clean continuum (lots of small sections)
clean_sections = [
    [4610.0, 4616.0], [4624.0, 4627.0], # between C II, N II, O II
    [4690.0, 4697.0], [4720.0, 4730.0], # between He I and Ar IV
    [4746.0, 4750.0], [4760.0, 4765.0], # between Fe III lines
    [4782.0, 4786.0], [4820.0, 4830.0], # next to Hb
    [4910.0, 4916.0], [5060.0, 5080.0], # to the red
    [5090.0, 5100.0], [5170.0, 5185.0], # to the red
]

cont_slices = []
for wavs in clean_sections:
    wavs = 1e-10*np.array(wavs)
    _, _, wpix = w.world_to_pixel_values([0, 0], [0, 0], wavs)
    cont_slices.append(slice(*wpix.astype(int)))


# Use median over each section to avoid weak lines
cont_maps = np.array([np.median(hdu.data[_, :, :], axis=0) for _ in cont_slices])
cont_wavpix = np.array([np.median(wavpix[_], axis=0) for _ in cont_slices])
# Inefficient but simple algorithm - loop over spaxels
bgdata = np.empty_like(hdu.data)
for j, i in itertools.product(range(ny), range(nx)):
    # Fit polynomial to BG
    try:
        p = T.fit(cont_wavpix, cont_maps[:, j, i], deg=2)
        # and fill in the BG spectrum of this spaxel
        bgdata[:, j, i] = p(wavpix)
    except:
        bgdata[:, j, i] = np.nan



for suffix, cube in [
        ["cont", bgdata],
        ["cont-sub", hdu.data - bgdata],
        # ["cont-div", hdu.data/bgdata],
]:
    outfile = infile.replace(".fits", f"-{suffix}.fits")
    fits.PrimaryHDU(header=hdu.header, data=cube).writeto(
        f"{OUTDIR}/{outfile}", overwrite=True)
    print(f"Written {outfile}")
