import sys
from pathlib import Path
import numpy as np
import pyneb
from astropy.io import fits

# Set up Blagrave 2007 extinction law
REDCORR = pyneb.extinction.red_corr.RedCorr(
    law='CCM89 Bal07', R_V=5.5, cHbeta=1.0)

def flambda(wav):
    """Find [(A_lam / A_Hb) - 1] as function of wavelength `wav`

    This is the same as given in Table 2 of Blagrave et al (2007)

    """
    return np.log10(REDCORR.getCorrHb(wav))


def CHb_from_RHbHa(RHbHa, balmer0=2.874):
    """Find base-10 extinction at H beta from balmer decrement `RHbHa`

    Assumes that the intrinsic Balmer decrement is `balmer0`

    """
    return np.log10(balmer0*RHbHa) / flambda(6563)


def CHb_from_R6563_9229(RBaPa, RBaPa0=112.0):
    """Find base-10 extinction at H beta from 6563/9229 decrement `RBaPa`

    Assumes that the intrinsic decrement is `RBaPa0`

    """
    return np.log10(RBaPa/RBaPa0) / (flambda(9229) - flambda(6563))



DATADIR = Path("../Orion-HH-data/MUSE")

if __name__ == '__main__':

    try:
        prefix = sys.argv[1]
        wav = int(sys.argv[2])
    except IndexError:
        print(f'Usage: {sys.argv[0]} LINEID WAV')

    hb_ha = fits.open(DATADIR / 'ratio-4861-6563.fits')[0].data
    chb = CHb_from_RHbHa(hb_ha)
    clam = (1.0 + flambda(wav))*chb
    fn = f"{prefix}.fits"
    hdu = fits.open(DATADIR / fn)[0]
    fn_new = f"{prefix}-excorr.fits"
    fits.PrimaryHDU(
        data=hdu.data*10**clam,
        header=hdu.header
    ).writeto(DATADIR / fn_new, overwrite=True)
    print(fn_new)
