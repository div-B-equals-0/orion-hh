TAB=[[38, 44, 4627.3, 4632.4, "", "N II 4631, [Ni II] 4628", "N_II-4631"], [44, 49, 4632.4, 4636.65, "", "N III 4634", "N_III-4634"], [49, 62, 4636.65, 4646.0, "4639, 4642", "N III 4641, N II 4643", "O_II-V1-4639-42"], [62, 71, 4646.0, 4654.5, "4649, 4651", "", "O_II-V1-4649-51"], [79, 83, 4662.15, 4664.7, 4662, "(Red wing of 4658)", "O_II-V1-4662"], [90, 102, 4671.5, 4680.85, "4674, 4676", "", "O_II-V1-4674-76"], [101, 119, 4679.15, 4694.45, "", "He II 4686 abs", "He_II-4686"]]
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
infile_bg = f"muse-hr-data-wavsec0{SUFFIX}-cont.fits"
hdu = fits.open(f"{DATADIR}/{infile}")["DATA"]
hdubg = fits.open(f"{DATADIR}/{infile_bg}")["DATA"]
w = WCS(hdu)
nwav, ny, nx = hdu.data.shape
wavpix = np.arange(nwav)
_, _, wavs = w.pixel_to_world_values([0]*nwav, [0]*nwav, wavpix)
wavs *= 1.e10
for i1, i2, _, _, _, _, label in TAB:
    outfile = f"{label}-sum{SUFFIX}.fits"
    imsum = np.sum(hdu.data[i1:i2, :, :], axis=0)
    fits.PrimaryHDU(
        header=w.celestial.to_header(),
        data=imsum
    ).writeto(f"{OUTDIR}/{outfile}", overwrite=True)
    print(f"Written {outfile}")

    outfile = f"{label}-wav{SUFFIX}.fits"
    imwav = np.sum(hdu.data[i1:i2, :, :]*wavs[i1:i2, None, None], axis=0)/imsum
    fits.PrimaryHDU(
        header=w.celestial.to_header(),
        data=imwav
    ).writeto(f"{OUTDIR}/{outfile}", overwrite=True)
    print(f"Written {outfile}")

    outfile = f"{label}-bg{SUFFIX}.fits"
    imsum = np.sum(hdubg.data[i1:i2, :, :], axis=0)
    fits.PrimaryHDU(
        header=w.celestial.to_header(),
        data=imsum
    ).writeto(f"{OUTDIR}/{outfile}", overwrite=True)
    print(f"Written {outfile}")
