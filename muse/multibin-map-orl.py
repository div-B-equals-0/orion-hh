import sys
from pathlib import Path
from distutils.dep_util import newer, newer_group
import numpy as np
sys.path.append("/Users/will/Dropbox/multibin-maps")
from rebin_utils import downsample, oversample, pad_array
from astropy.io import fits

nlist = [1, 2, 4, 8, 16, 32, 64, 128, 256]
mingoods = [2, 2, 2, 2, 2, 2, 2, 2, 2]


try:
    datadir = Path(sys.argv[1])
    infile = sys.argv[2]
    contfile = sys.argv[3]
    cont_threshold = float(sys.argv[4])
except:
    sys.exit('Usage: {} DATADIR INFILE CONTFILE CONT_THRESHOLD'.format(sys.argv[0]))


hdu = fits.open(datadir / infile)[0]
if hdu.data is None:
    hdu = fits.open(datadir / infile)[1]
hdr = hdu.header
# Maximum binning
nmax = nlist[-1]

# continuum image
chdu = fits.open(datadir / contfile)[0]
if chdu.data is None:
    chdu = fits.open(datadir / contfile)[1]

# Pad arrays to nearest multiple of nmax
im = pad_array(hdu.data, nmax)
cim = pad_array(chdu.data, nmax)

w = np.ones_like(im)
starmask = cim > cont_threshold

# If we pad the starmask and combine it with the padded image, then we
# automatically deal with the case where the input files have already
# been padded
m =  np.isfinite(im) & (~pad_array(starmask, nmax))

for n, mingood in zip(nlist, mingoods):
    im[~m] = 0.0
    outfile = infile.replace('.fits', '-bin{:03d}.fits'.format(n))
    if n == nlist[0]:
        # Do dependency checking on the first iteration
        if not newer(datadir / infile, datadir / outfile):
            # Bail out if dependency not newer than target
            sys.exit(outfile + ' is already up to date.')
    print('Saving', outfile)
    # Save both the scaled image and the weights, but at the full resolution
    fits.HDUList([
        fits.PrimaryHDU(),
        fits.ImageHDU(data=oversample(im, n), header=hdr, name='scaled'),
        fits.ImageHDU(data=oversample(w, n), header=hdr, name='weight'),
    ]).writeto(datadir / outfile, clobber=True)
    # Now do the rebinning by a factor of two
    [im,], m, w = downsample([im,], m, weights=w, mingood=mingood)
