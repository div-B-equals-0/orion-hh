# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.7
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %pylab inline

from astropy.io import fits
import numpy as np
from wfc3_utils import get_filter, get_interpolated_filter
from astropy.table import Table, pprint
pprint.MAX_WIDTH.set(120)

# The `table_M42_sr.fits` file has the metadata for each pixel - position

thdu = fits.open('Manu-Data/M42_red/table_M42_lr.fits')[1]

Table(thdu.data)

m = np.hypot(thdu.data.dRA - 15.0 , thdu.data.dDEC + 70.0) < 10.0

# This is a new mask, centered on the O II peak in the Big Arc

sum(m)

Table(thdu.data[m])

hdu = fits.open('Manu-Data/M42_red/M42_lr.fits')[0]
shdu = fits.open('Manu-Data/M42_red/M42_sr.fits')[0]

hdu.header

hdu.data.shape, m.shape

spectra = hdu.data[m]*thdu.data[m].factor2[:, None]
spectra_s = shdu.data[m]*thdu.data[m].factor2[:, None]

spectra.shape

spectra[0]

wav0, dwav, nwav = (hdu.header[k] for k in ('CRVAL1', 'CDELT1', 'NAXIS1'))

wavs = wav0 + np.arange(nwav)*dwav

wavs

# +
wavs_656, f656n = get_filter("F656N", return_wavelength=True)
wavs_658, f658n = get_filter("F658N", return_wavelength=True)

plot(wavs, (wavs/6000)*spectra.mean(axis=0), 'k', label="long")
plot(wavs, (wavs/6000)*spectra_s.mean(axis=0), 'r', lw=0.3, label="short")
for spec in spectra:
    plot(wavs, (wavs/6000)*spec, 'k', alpha=0.1)
ylim(0.0, 2e-14)
xticks(np.arange(5600.0, 7000.0, 100.0))
grid()
legend()
xlabel("Wavelength, Angstrom")
ylabel("lambda F_lambda, erg/cm2/s")
#yscale('log')
gcf().set_size_inches(12,8)
gcf().savefig("manuel-red.pdf")

# +
spec_combine = np.zeros_like(spectra)
shortmask = (np.abs(wavs - 6563.0) < 5.5) | (np.abs(wavs - 6584.0) < 3.0) 
spec_combine = spectra_s
#spec_combine[shortmask[None,:]] = spectra_s[shortmask[None,:]]

plot(wavs, spec_combine.mean(axis=0), 'k', label="long")
for spec in spec_combine:
    plot(wavs, spec, 'k', alpha=0.3)
    
for spec in spectra:
    plot(wavs, np.where(shortmask, np.nan, spec), 'r', alpha=0.3)

plot(wavs_658, 1e-11*f658n, 'k--', lw=5, alpha=0.2)
plot(wavs_656, 1e-11*f656n, 'k--', lw=5, alpha=0.2)
#yscale('log')
ylim(1.0e-15, 1.0e-11)
xlim(6530.0, 6620.0)
#xticks(np.arange(5600.0, 6900.0, 100.0))
grid()
yscale('log')
legend()
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
gcf().set_size_inches(12,8)
# -

plot(wavs, spectra.mean(axis=0), 'k', label="long")
plot(wavs, spectra_s.mean(axis=0), 'r', label="short")
#yscale('log')
ylim(0.0, 0.15e-13)
xlim(5700.0, 5800.0)
#xticks(np.arange(5600.0, 6900.0, 100.0))
grid()
legend()
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
gcf().set_size_inches(12,8)

# ## Now look at the green and blue spectra too

thdu_g = fits.open('Manu-Data/M42_green/table_M42_lg.fits')[1]
thdu_b = fits.open('Manu-Data/M42_blue/table_M42_lb.fits')[1]

# +
mb = hypot(thdu_b.data.dRA - 15.0 , thdu_b.data.dDEC + 70.0) < 10.0
mg = hypot(thdu_g.data.dRA - 15.0 , thdu_g.data.dDEC + 70.0) < 10.0

Table(thdu_b.data[mb])
# -

Table(thdu_g.data[mg])

# Check that there is no difference between the fibers used in the long and short spectra

# +
thdu_sg = fits.open('Manu-Data/M42_green/table_M42_sg.fits')[1]
thdu_sb = fits.open('Manu-Data/M42_blue/table_M42_sb.fits')[1]
thdu_sr = fits.open('Manu-Data/M42_red/table_M42_sr.fits')[1]
msb = hypot(thdu_sb.data.dRA - 15.0 , thdu_sb.data.dDEC + 70.0) < 10.0
msg = hypot(thdu_sg.data.dRA - 15.0 , thdu_sg.data.dDEC + 70.0) < 10.0

assert( np.alltrue(thdu_sb.data[msb].id_ap == thdu_b.data[mb].id_ap) )
assert( np.alltrue(thdu_sg.data[msg].id_ap == thdu_g.data[mg].id_ap) )
# -

# Yep, they are all the same thankfully

hdu_b = fits.open('Manu-Data/M42_blue/M42_lb.fits')[0]
shdu_b = fits.open('Manu-Data/M42_blue/M42_sb.fits')[0]
hdu_g = fits.open('Manu-Data/M42_green/M42_lg.fits')[0]
shdu_g = fits.open('Manu-Data/M42_green/M42_sg.fits')[0]

spectra_b = hdu_b.data[mb]*thdu_b.data[mb].factor2[:, None]
spectra_sb = shdu_b.data[msb]*thdu_sb.data[msb].factor2[:, None]
spectra_g = hdu_g.data[mg]*thdu_g.data[mg].factor2[:, None]
spectra_sg = shdu_g.data[msg]*thdu_sg.data[msg].factor2[:, None]

wav0, dwav, nwav = (hdu_b.header[k] for k in ('CRVAL1', 'CDELT1', 'NAXIS1'))
wavs_b = wav0 + np.arange(nwav)*dwav
wav0, dwav, nwav = (hdu_g.header[k] for k in ('CRVAL1', 'CDELT1', 'NAXIS1'))
wavs_g = wav0 + np.arange(nwav)*dwav

# The green spectrum on its own.  I mask out the long spectrum for $ F > 10^{-13}$ since it goes funny when it gets saturated.
#

msat = spectra_g.mean(axis=0) < 1.e-13
plot(wavs_g[msat], spectra_g.mean(axis=0)[msat], 'k', label="long")
plot(wavs_g, spectra_sg.mean(axis=0), 'g', lw=0.3, label="short")
yscale('log')
ylim(4e-15, 5e-14)
#ylim(0.0, 2e-14)
#xticks(np.arange(5600.0, 6900.0, 100.0))
grid()
legend()
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
gcf().set_size_inches(12,8)

# The blue spectrum on its own

msat = spectra_b.mean(axis=0) < 1.e-11
plot(wavs_b[msat], spectra_b.mean(axis=0)[msat], 'k', label="long")
plot(wavs_b, spectra_sb.mean(axis=0), 'b', lw=0.3, label="short")
yscale('log')
ylim(4e-15, 7e-14)
#xticks(np.arange(5600.0, 6900.0, 100.0))
grid()
legend()
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
gcf().set_size_inches(12,8)

# Zooms on the blue range.

# +
f469n = get_interpolated_filter("F469N", wavs_b)
plot(wavs_b, 1e-13*f469n, '--k', alpha=0.5, label='F469N')
fill_between(wavs_b, f469n*spectra_b.mean(axis=0)/f469n.max(), np.zeros_like(wavs_b), 
             alpha=0.25, zorder=-100, facecolor='k')

plot(wavs_b[msat], spectra_b.mean(axis=0)[msat], 'k', label="long")
plot(wavs_b, spectra_sb.mean(axis=0), 'b', lw=0.3, label="short")
plot(wavs_g, 0.91*spectra_g.mean(axis=0), 'g', label="long green")
#yscale('log')
xlim(4450, 4850)
ylim(0.0, 0.3e-13)
minorticks_on()
grid(ls='-', c='b', lw=0.6, alpha=0.3)
grid(ls='-', c='b', lw=0.6, alpha=0.05, which='minor')
legend()
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
gcf().set_size_inches(12,8)
tight_layout()
# -

for spec, id_ in zip(spectra_b, thdu_b.data[mb].id_ap):
    cont_mask = (wavs_b>=4580) & (wavs_b<=4590)
    norm = spec[cont_mask].mean()
    plot(wavs_b, spec/norm, lw=0.6, label="fiber {}".format(id_), alpha=0.1)
#yscale('log')
mspec = spectra_b.mean(axis=0)
mspec /= mspec[cont_mask].mean()
plot(wavs_b, mspec, "k", lw=2)
xlim(4580, 4740)
ylim(0.8, 1.2)
minorticks_on()
grid(ls='-', c='b', lw=0.6, alpha=0.3)
grid(ls='-', c='b', lw=0.6, alpha=0.05, which='minor')
#legend(ncol=5, loc="lower left")
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
gcf().set_size_inches(12,8)

# Here we have absorption in the 4686 line.  

for spec, id_ in zip(spectra_b[::5], thdu_b.data[mb].id_ap[::5]):
    cont_mask = (wavs_b>=4520) & (wavs_b<=4530)
    norm = spec[cont_mask].mean()
    plot(wavs_b, spec/norm, lw=0.6, label="fiber {}".format(id_))
#yscale('log')
mspec = spectra_b.mean(axis=0)
mspec /= mspec[cont_mask].mean()
plot(wavs_b, mspec, "k", lw=2)
xlim(4480, 4580)
ylim(0.8, 1.2)
minorticks_on()
grid(ls='-', c='b', lw=0.6, alpha=0.3)
grid(ls='-', c='b', lw=0.6, alpha=0.05, which='minor')
legend(ncol=5, loc="lower left")
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
gcf().set_size_inches(12,8)

# Here we have absorption in the He II 4541 line at a level of about 4%.  The stellar spectrum of th1C shows absorption at the 20% level. 

for spec, id_ in zip(spectra_b, thdu_b.data[mb].id_ap):
    cont_mask = (wavs_b>=4395) & (wavs_b<=4405)
    norm = spec[cont_mask].mean()
    plot(wavs_b, spec/norm, lw=0.6, label="fiber {}".format(id_), alpha=0.1)
#yscale('log')
mspec = spectra_b.mean(axis=0)
mspec /= mspec[cont_mask].mean()
plot(wavs_b, mspec, "k", lw=2)
xlim(4300, 4600)
ylim(0.7, 1.2)
minorticks_on()
grid(ls='-', c='b', lw=0.6, alpha=0.3)
grid(ls='-', c='b', lw=0.6, alpha=0.05, which='minor')
#legend(ncol=5, loc="lower left")
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
gcf().set_size_inches(12,8)

# There will be absorption underneath the two He I lines 4387 and 4471, but it is impossible to see because the emission swamps it. 

for spec, id_ in zip(spectra_b, thdu_b.data[mb].id_ap):
    cont_mask = (wavs_b>=4220) & (wavs_b<=4240)
    norm = spec[cont_mask].mean()
    plot(wavs_b, spec/norm, lw=0.6, label="fiber {}".format(id_), alpha=0.1)
#yscale('log')
mspec = spectra_b.mean(axis=0)
mspec /= mspec[cont_mask].mean()
plot(wavs_b, mspec, "k", lw=2)
xlim(4100, 4300)
ylim(0.6, 1.6)
minorticks_on()
grid(ls='-', c='b', lw=0.6, alpha=0.3)
grid(ls='-', c='b', lw=0.6, alpha=0.05, which='minor')
#legend(ncol=5, loc="lower left")
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
title("Zoom in on region around O II V36 4189")
gcf().set_size_inches(12,8)

# We can see the stellar absorption at He II 4199.83, with a depth of about 2-3%.  From Simon-Diaz et al (2006) the absorption from th1C is about 12%, implying that 25% of the continuum is scattered th1C starlight. 
#   

for spec, id_ in zip(spectra_b, thdu_b.data[mb].id_ap):
    cont_mask = (wavs_b>=4050) & (wavs_b<=4060)
    norm = spec[cont_mask].mean()
    plot(wavs_b, spec/norm, lw=0.6, label="fiber {}".format(id_), alpha=0.1)
#yscale('log')
mspec = spectra_b.mean(axis=0)
mspec /= mspec[cont_mask].mean()
plot(wavs_b, mspec, "k", lw=2)
xlim(4000, 4140)
ylim(0.6, 1.8)
minorticks_on()
grid(ls='-', c='b', lw=0.6, alpha=0.3)
grid(ls='-', c='b', lw=0.6, alpha=0.05, which='minor')
#legend(ncol=5, loc="lower left")
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
title("Zoom in on region around O II V48a 4089")
gcf().set_size_inches(12,8)

for spec, id_ in zip(spectra_g, thdu_g.data[mg].id_ap):
    cont_mask = (wavs_g>=4500) & (wavs_g<=4530)
    norm = spec[cont_mask].mean()
    plot(wavs_g, spec/norm, lw=0.6, label="fiber {}".format(id_), alpha=0.1)
#yscale('log')
mspec = spectra_g.mean(axis=0)
mspec /= mspec[cont_mask].mean()
plot(wavs_g, mspec, "k", lw=2)
xlim(4500, 4700)
ylim(0.7, 1.2)
minorticks_on()
grid(ls='-', c='b', lw=0.6, alpha=0.3)
grid(ls='-', c='b', lw=0.6, alpha=0.05, which='minor')
#legend(ncol=5, loc="lower left")
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
gcf().set_size_inches(12,8)

for spec, id_ in zip(spectra_b, thdu_b.data[mb].id_ap):
    cont_mask = (wavs_b>=3980) & (wavs_b<=4000)
    norm = spec[cont_mask].mean()
    plot(wavs_b, spec/norm, lw=0.6, label="fiber {}".format(id_), alpha=0.1)
#yscale('log')
mspec = np.median(spectra_b, axis=0)
mspec /= mspec[cont_mask].mean()
plot(wavs_b, mspec, "k", lw=2)
plot(wavs_b, 0.8 + (mspec-0.8)/30, "r", lw=1, alpha=0.6)
xlim(3800, 4000)
ylim(0.7, 1.5)
minorticks_on()
grid(ls='-', c='b', lw=0.6, alpha=0.3)
grid(ls='-', c='b', lw=0.6, alpha=0.1, which='minor')
#legend(ncol=5, loc="lower left")
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
title("[Ne III] 3869 and and He I 3889 (also Si II 3863 and He I 3878)")
gcf().set_size_inches(12,8)

# [S III] 3798 is hidden by the H line, so this is useless

# +
for spec, id_ in zip(spectra_b, thdu_b.data[mb].id_ap):
    cont_mask = (wavs_b>=3808) & (wavs_b<=3815)
    norm = spec[cont_mask].mean()
    plot(wavs_b, spec/norm, lw=0.6, label="fiber {}".format(id_), alpha=0.1)

mspec = np.mean(spectra_b, axis=0)
mspec /= mspec[cont_mask].mean()
plot(wavs_b, mspec, "k", lw=2)
plot(wavs_b, 0.8 + (mspec-0.8)/10, "r", lw=1, alpha=0.6)
#yscale('log')
xlim(3750, 3850)
ylim(0.5, 1.5)
minorticks_on()
grid(ls='-', c='b', lw=0.6, alpha=0.3)
grid(ls='-', c='b', lw=0.6, alpha=0.05, which='minor')
#legend(ncol=4, loc="upper left")
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
title("[S III]+H10 3798 in middle, H11 3771 to blue, He I 3820 and H9 3835 to red")
gcf().set_size_inches(12,8)
# -

for spec, id_ in zip(spectra_b, thdu_b.data[mb].id_ap):
    cont_mask = (wavs_b>=3640) & (wavs_b<=3660)
    norm = spec[cont_mask].mean()
    plot(wavs_b, spec/norm, lw=0.6, label="fiber {}".format(id_), alpha=0.1)
#yscale('log')
mspec = np.mean(spectra_b, axis=0)
mspec /= mspec[cont_mask].mean()
plot(wavs_b, mspec, "k", lw=2)
plot(wavs_b, 0.4 + (mspec-0.4)/30, "r", lw=1, alpha=0.6)
xlim(3560, 3760)
ylim(0.3, 1.5)
minorticks_on()
grid(ls='-', c='b', lw=0.6, alpha=0.3)
grid(ls='-', c='b', lw=0.6, alpha=0.05, which='minor')
#legend(ncol=4, loc="upper left")
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
title("[O II] blue lines and Balmer limit")
gcf().set_size_inches(12,8)

for spec, id_ in zip(spectra_b, thdu_b.data[mb].id_ap):
    cont_mask = (wavs_b>=3640) & (wavs_b<=3660)
    norm = spec[cont_mask].mean()
    plot(wavs_b, spec/norm, lw=0.6, label="fiber {}".format(id_), alpha=0.1)
#yscale('log')
mspec = np.mean(spectra_b, axis=0)
mspec /= mspec[cont_mask].mean()
plot(wavs_b, mspec, "k", lw=2)
xlim(3500, 4200)
ylim(0.0, 2.5)
minorticks_on()
grid(ls='-', c='b', lw=0.6, alpha=0.3)
grid(ls='-', c='b', lw=0.6, alpha=0.05, which='minor')
#legend(ncol=4, loc="upper left")
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
title("Zoom out on the Balmer jump")
gcf().set_size_inches(12,8)

for spec, id_ in zip(spectra_b, thdu_b.data[mb].id_ap):
    cont_mask = (wavs_b>=4720) & (wavs_b<=4730)
    norm = spec[cont_mask].mean()
    plot(wavs_b, spec/norm, lw=0.6, label="fiber {}".format(id_), alpha=0.1)
#yscale('log')
mspec = np.median(spectra_b, axis=0)
mspec /= mspec[cont_mask].mean()
plot(wavs_b, mspec, "k", lw=2)
plot(wavs_b, 0.4 + (mspec-0.4)/30, "r", lw=2, alpha=0.8)
xlim(4720, 4950)
ylim(0.0, 1.3)
minorticks_on()
grid(ls='-', c='b', lw=0.6, alpha=0.3)
grid(ls='-', c='b', lw=0.6, alpha=0.05, which='minor')
#legend(ncol=4, loc="upper left")
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
title("Going back to the 4700 end")
gcf().set_size_inches(12,8)

for spec, id_ in zip(spectra_b, thdu_b.data[mb].id_ap):
    cont_mask = (wavs_b>=3760) & (wavs_b<=3840)
    norm = spec[cont_mask].mean()
    plot(wavs_b, spec/norm, lw=0.6, label="fiber {}".format(id_))
#yscale('log')
xlim(3900, 3960)
ylim(0, 2)
minorticks_on()
grid(ls='-', c='b', lw=0.6, alpha=0.3)
grid(ls='-', c='b', lw=0.6, alpha=0.05, which='minor')
legend(ncol=4, loc="upper left")
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
title("C II 3920 doublet and He I 3927")
gcf().set_size_inches(12,8)

# ### Red–green overlap

wavs_547, f547m = get_filter("F547M", return_wavelength=True)
wavs_575, fq575n = get_filter("FQ575N", return_wavelength=True)

# +
msat_r = spectra.mean(axis=0) < 1.e-13
msat_g = spectra_g.mean(axis=0) < 1.e-13
plot(wavs_547, 1e-13*f547m, '--k', alpha=0.5, label='F547M')
fill_between(wavs_547, 1e-13*f547m, np.zeros_like(wavs_547), alpha=0.05, zorder=-100, facecolor='k')
plot(wavs_575, 1e-13*fq575n, ':k', alpha=0.8, label='FQ575N')
fill_between(wavs_575, 1e-13*fq575n, np.zeros_like(wavs_575), alpha=0.1, zorder=-50, facecolor='k')
plot(wavs[msat_r], spectra.mean(axis=0)[msat_r], 'k', label="long red")
plot(wavs, spectra_s.mean(axis=0), 'r', lw=0.3, label="short red")
plot(wavs_g[msat_g], spectra_g.mean(axis=0)[msat_g], 'k', label="long green")
plot(wavs_g, spectra_sg.mean(axis=0), 'g', lw=0.3, label="short green")
#yscale('log')
xlim(5000.0, 6000.0)
ylim(0.0, 0.3e-13)
xticks(np.arange(5000.0, 6100.0, 100.0))
minorticks_on()
grid(ls='-', c='b', lw=0.3, alpha=0.3)
grid(ls='-', c='b', lw=0.3, alpha=0.05, which='minor')

#grid(ls='-', c='b', lw=0.5)
legend()
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
title("Red/green overlap region for mean spectrum within 5'' of (-30'', -30'')")
gcf().set_size_inches(12,4)
gcf().savefig("manu-red-green-comparison.pdf")
# -

# Compare the positions of the fibers:

scatter(thdu.data[m].dRA, thdu.data[m].dDEC, s=300.0, alpha=0.5, facecolor='r')
scatter(thdu_g.data[mg].dRA, thdu_g.data[mg].dDEC, s=300.0, alpha=0.5, facecolor='g')
scatter(thdu_b.data[mb].dRA, thdu_b.data[mb].dDEC, s=300.0, alpha=0.5, facecolor='b')
xlim(0.0, -40.0)
ylim(-40.0, 0.0)
xlabel("dRA, arcsec")
ylabel("dDEC, arcsec")
title("All fibers within 5'' of (-30'', -30'')")
grid()
gcf().set_size_inches(6,6)
gcf().savefig("manu-fiber-positions.pdf")

# Again but with $\lambda F_\lambda$

# +
msat_r = spectra.mean(axis=0) < 1.e-13
msat_g = spectra_g.mean(axis=0) < 1.e-13
plot(wavs_547, 5500*1e-13*f547m, '--k', alpha=0.5, label='F547M')
fill_between(wavs_547, 5500*1e-13*f547m, np.zeros_like(wavs_547), alpha=0.05, zorder=-100, facecolor='k')
plot(wavs_575, 5500*1e-13*fq575n, ':k', alpha=0.8, label='FQ575N')
fill_between(wavs_575, 5500*1e-13*fq575n, np.zeros_like(wavs_575), alpha=0.1, zorder=-50, facecolor='k')
plot(wavs[msat_r], wavs[msat_r]*spectra.mean(axis=0)[msat_r], 'k', label="long red")
plot(wavs, wavs*spectra_s.mean(axis=0), 'r', lw=0.3, label="short red")
plot(wavs_g[msat_g], wavs_g[msat_g]*spectra_g.mean(axis=0)[msat_g], 'k', label="long green")
plot(wavs_g, wavs_g*spectra_sg.mean(axis=0), 'g', lw=0.3, label="short green")
plot(wavs_547, 5e-11*np.ones_like(wavs_547), ':k', lw=2, alpha=0.5)
#yscale('log')
xlim(5000.0, 6000.0)
ylim(0.0, 5500*0.3e-13)
xticks(np.arange(5000.0, 6100.0, 100.0))
minorticks_on()
grid(ls='-', c='b', lw=0.3, alpha=0.3)
grid(ls='-', c='b', lw=0.3, alpha=0.05, which='minor')

#grid(ls='-', c='b', lw=0.5)
legend()
xlabel("Wavelength, Angstrom")
ylabel("lambda F_lambda, erg/cm2/s")
title("Red/green overlap region for mean spectrum within 5'' of (-30'', -30'')")
gcf().set_size_inches(12,4)
gcf().savefig("manu-red-green-lamFlam-comparison.pdf")
# -

# ### Choosing the best fibers

# Manu says the following: 

# *Por último, por donde están los peores ajustes tras la reducción son por el lado al rojo de cada rango espectral, así por ejemplo en el rango verde yo cortaría en 5600. Aunque en tu gráfico parece que tus fibras no están muy afectadas de tal efecto, las fibras más afectadas tienen un id_ap entre 1-61, si no recuerdo mal (prueba a imprimir esa columna usando la máscara a ver que valores tienes, si están entre 61 y 300 son las mejores fibras).*



# +
fig, (ax1, ax2) = subplots(1, 2, sharey=True)
size = 15.0
alpha = 0.25


badmask = thdu_b.data.id_ap <= 61
ax1.scatter(thdu_b.data[badmask].dRA, thdu_b.data[badmask].dDEC, 
            s=size, alpha=alpha, facecolor='b', linewidths=0.0)
ax2.scatter(thdu_b.data[~badmask].dRA, thdu_b.data[~badmask].dDEC, 
            s=size, alpha=alpha, facecolor='b', linewidths=0.0)

badmask = thdu_g.data.id_ap <= 61
ax1.scatter(thdu_g.data[badmask].dRA, thdu_g.data[badmask].dDEC, 
            s=size, alpha=alpha, facecolor='g', linewidths=0.0)
ax2.scatter(thdu_g.data[~badmask].dRA, thdu_g.data[~badmask].dDEC, 
            s=size, alpha=alpha, facecolor='g', linewidths=0.0)

badmask = thdu.data.id_ap <= 61
ax1.scatter(thdu.data[badmask].dRA, thdu.data[badmask].dDEC, 
            s=size, alpha=alpha, facecolor='r', linewidths=0.0)
ax2.scatter(thdu.data[~badmask].dRA, thdu.data[~badmask].dDEC, 
            s=size, alpha=alpha, facecolor='r', linewidths=0.0)

ax1.set_xlim(-150.0, 150.0)
ax1.set_ylim(-200.0, 200.0)
ax1.set_xlabel("dRA, arcsec")
ax1.set_ylabel("dDEC, arcsec")
ax1.set_title("Fibers with calibration problems at red end")
ax1.grid()

ax2.set_xlim(-150.0, 150.0)
ax2.set_ylim(-200.0, 200.0)
ax2.set_title("Fibers with the best calibration")
ax2.grid()
ax2.set_xlabel("dRA, arcsec")

fig.set_size_inches(12,8)
fig.tight_layout()
fig.savefig("manu-all-fiber-positions.pdf")
fig.savefig("manu-all-fiber-positions.jpg", dpi=300)
# -

# ### Big sum of everything

spectra_b = hdu_b.data*thdu_b.data.factor2[:, None]
spectra_sb = shdu_b.data*thdu_sb.data.factor2[:, None]*thdu_sb.data.factor1[:, None]
spectra_g = hdu_g.data*thdu_g.data.factor2[:, None]
spectra_sg = shdu_g.data*thdu_sg.data.factor2[:, None]*thdu_sg.data.factor1[:, None]

spectra_r = hdu.data*thdu.data.factor2[:, None]
spectra_sr = shdu.data*thdu_sr.data.factor2[:, None]*thdu_sr.data.factor1[:, None]

radius = 15.0
msr = hypot(thdu_sr.data.dRA + 30.0 , thdu_sr.data.dDEC + 30.0) < 15.0
msg = hypot(thdu_sg.data.dRA + 30.0 , thdu_sg.data.dDEC + 30.0) < 15.0
msb = hypot(thdu_sb.data.dRA + 30.0 , thdu_sb.data.dDEC + 30.0) < 15.0
plot(wavs, spectra_r[msr].mean(axis=0))
plot(wavs_g, spectra_g[msg].mean(axis=0))
plot(wavs, spectra_sr[msr].mean(axis=0))
plot(wavs_g, spectra_sg[msg].mean(axis=0))
ylim(0, 3e-14)
xlim(5500, 5950)
grid()
minorticks_on()
grid(which='minor', ls='-', alpha=0.1)
gcf().set_size_inches(12, 6)

plot(wavs_b, spectra_b.sum(axis=0))
ylim(0, 2e-10)
xlim(4500, 4800)
grid()
minorticks_on()
grid(which='minor', ls='-', alpha=0.1)
gcf().set_size_inches(12, 6)

ratio = spectra_sr/spectra_r
for r in ratio[msr]:
    plot(wavs, r, 'k', lw=0.1, alpha=0.1)
ylim(0.0, 2.0)
gcf().set_size_inches(12,6)

medrat = np.median(ratio, axis=1)
hist(medrat, bins=200, range=(0.0, 2.0))

ratio = spectra_sg/spectra_g
for r in ratio[msg]:
    plot(wavs_g, r, 'k', lw=0.1, alpha=0.1)
ylim(0.0, 2.0)
gcf().set_size_inches(12,6)

medrat = np.median(ratio, axis=1)
hist(medrat, bins=200, range=(0.0, 2.0))

ratio = spectra_sb/spectra_b
for r in ratio[msb]:
    plot(wavs_b, r, 'k', lw=0.1, alpha=0.2)
ylim(0.0, 2.0)
gcf().set_size_inches(12,6)

medrat = np.median(ratio, axis=1)
hist(medrat, bins=200, range=(0.0, 2.0))

# ## Checking out the calibration data

tab = Table.read("manu_compare_FQ575N_F547M.tab", format="ascii.tab", fill_values=('--', 0.0))

m = ~tab['E5755'].mask & ~tab['FQ575N'].mask & ~tab['F547M'].mask
scatter(tab['x'][m], tab['y'][m], s=250.0, c=tab['E5755'][m], 
        alpha=0.4, vmin=0, vmax=15)
colorbar()
xlim(0.0, -90.0)
ylim(-90.0, 0.0)rm
axis('equal')
gcf().set_size_inches(12,12)

scatter(tab['x'][m], tab['y'][m], c=tab['FQ575N'][m]/tab['F547M'][m], 
        s=250.0, alpha=0.4, vmin=0.025, vmax=0.045)
colorbar()
xlim(0.0, -90.0)
ylim(-90.0, 0.0)
axis('equal')
gcf().set_size_inches(12,12)

tab['E5755']

# ## The O II recombination line density

tabg = Table.read("manu-diag-green-fluxes.tab", format="ascii.tab")
tabb = Table.read("manu-diag-blue-fluxes.tab", format="ascii.tab")

# Correct the blue spectra in the lower right-hand corner.

mask = np.hypot(tabb['x'], tabb['y']) > 85.0
for k in tabb.colnames:
    if k.startswith('F'):
        tabb[k][mask] /= 1.3

tabg.colnames

tabb.colnames


def find_xytotrat(tab):
    x = tab['x']
    y = tab['y']
    total = tab['F4649'] + tab['F4639'] + tab['F4651'] + tab['F4662']
    ratio = tab['F4649']/(tab['F4639'] + tab['F4651'] + tab['F4662'])
    return x, y, total, ratio


xb, yb, totalb, ratiob = find_xytotrat(tabb)
xg, yg, totalg, ratiog = find_xytotrat(tabg)

# Histogram of O II density ratio, from the blue spectra

hist(ratiob, range=[0.0, 1.2])

# Histogram of O II density ratio, from the green spectra

hist(ratiog, range=[0.0, 1.2])

hist(totalb, range=[0.0, 30.0])

hist(totalg, range=[0.0, 30.0])

# Maps of O II 4650 density ratio from blue spectra (left) and green spectra (right).  The map looks a lot better from the blue spectra. 

fig, (ax1, ax2) = subplots(1, 2, sharex=True, sharey=True)
scat1 = ax1.scatter(xb, yb, c=ratiob, vmin=0.4, vmax=1.1, s=200, cmap=cm.jet, alpha=0.4)
scat2 = ax2.scatter(xg, yg, c=ratiog, vmin=0.4, vmax=1.1, s=200, cmap=cm.jet, alpha=0.4)
fig.colorbar(scat1, ax=ax1)
fig.colorbar(scat2, ax=ax2)
ax1.set_xlim(0.0, -90.0)
ax1.set_ylim(-90.0, 0.0)
ax1.axis('equal')
ax2.axis('equal')
fig.set_size_inches(12,6)
fig.tight_layout()

fig, (ax1, ax2) = subplots(1, 2, sharex=True, sharey=True)
scat1 = ax1.scatter(xb, yb, c=totalb, vmin=0.0, vmax=20.0, s=200, cmap=cm.jet, alpha=0.8)
scat2 = ax2.scatter(xg, yg, c=totalg, vmin=0.0, vmax=20.0, s=200, cmap=cm.jet, alpha=0.8)
fig.colorbar(scat1, ax=ax1)
fig.colorbar(scat2, ax=ax2)
ax1.set_xlim(0.0, -90.0)
ax1.set_ylim(-90.0, 0.0)
ax1.axis('equal')
ax2.axis('equal')
fig.set_size_inches(12,6)
fig.tight_layout()

clrat = tabg["F5518"]/tabg["F5538"]
arrat = tabg["F4711"]/tabg["F4740"]

hist(clrat)

hist(1./arrat, range=[0.0, 1.])

fig, (ax1, ax2) = subplots(1, 2, sharex=True, sharey=True)
scat1 = ax1.scatter(xg, yg, c=clrat, vmin=0.4, vmax=1.1, s=200, cmap=cm.jet_r, alpha=0.8)
scat2 = ax2.scatter(xg, yg, c=1./arrat, vmin=0.0, vmax=1.2, s=200, cmap=cm.jet, alpha=0.8)
fig.colorbar(scat1, ax=ax1)
fig.colorbar(scat2, ax=ax2)
ax1.set_xlim(0.0, -90.0)
ax1.set_ylim(-90.0, 0.0)
ax1.axis('equal')
ax1.set_title('[Cl III] 5518/5538 green')
ax2.set_title('[Ar IV] 4711/4740 green')
ax2.axis('equal')
fig.set_size_inches(12,6)
fig.tight_layout()

# Plot of O II ratio versus [Cl III] ratio for the green spectra.

scatter(clrat, ratiog, alpha=0.5)
xlim(0.4, 1.1)
xlabel("[Cl III] ratio")
ylabel("O II ratio")
ylim(0.3, 1.3)

# There is some correlation there, but it is very noisy.

# Try to do the same for the O II ratio from the blue spectra.  This requires interpolation. 

from scipy.interpolate import griddata

xi = np.vstack([xb, yb]).T

xi.shape

ratio_bg = griddata((xb, yb), ratiob, (xg, yg), method='linear')
total_bg = griddata((xb, yb), totalb, (xg, yg), method='linear')


# +
def oii_den_ratio(den):
    x = np.log10(den)
    return 0.65*(np.tanh(1.2*(x - 3.45)) + 1.2)

def cliii_den_ratio(den):
    x = np.log10(den)
    return 0.55*(np.tanh(1.2*(3.55 - x)) + 1.62)

denrange = np.logspace(2.0, 5.0, 200)

# -

scatter(clrat, ratio_bg, c=radius_g, 
        vmin=30, vmax=100, s=100, alpha=0.6)
plot(cliii_den_ratio(denrange), oii_den_ratio(denrange), 
     '-k', label='n([Cl III]) = n(O II)')
plot(cliii_den_ratio(2*denrange), oii_den_ratio(denrange), 
     '--k', label='n([Cl III]) = 2 x n(O II)')
plot(cliii_den_ratio(4*denrange), oii_den_ratio(denrange), 
     ':k', label='n([Cl III]) = 4 x n(O II)')
xlim(0.4, 1.1)
xlabel("[Cl III] ratio: 5518 / 5538")
ylabel("O II ratio: 4649 / (4639 + 4651 + 4662)")
ylim(0.4, 1.1)
cb = colorbar()
cb.set_label('Radius, arcsec')
legend()
gcf().set_size_inches(9, 8)
savefig("oii-vs-cliii-densities.pdf")

# ### Checking ratios of recomb lines that should be constant and/or T- or N-sensitive

# The O II 4642/(4639+49+51+62) ratio should be roughly constant at 0.32, except for the fact that 4642 may be blended with N III

hist( tabb['F4642']/totalb, bins=100, range=[0.0, 2.0])
xlabel("O II 4642 / (4639 + 4649 + 4651 + 4662)")
ylabel("Number of apretures")

radius_b = np.hypot(xb, yb)
radius_g = np.hypot(xg, yg)

scatter(totalb, tabb['F4642'], c="r", alpha=0.8)
scatter(totalb, tabb['F4676'], c="b", alpha=0.8)
xlim(0.0, 30.0)
ylim(0.0, 10.0)
xgrid = np.linspace(0.0, 30.0)
plot(xgrid, 0.32*xgrid, "r--", label="O II 4642 theoretical slope: 0.32")
plot(xgrid, 0.32*1.43*xgrid, "r-", label="O II 4642 plus 43% contamination")
plot(xgrid, 0.15*xgrid, "b-", label="O II (4674 + 4676) theoretical slope: 0.15")
xlabel("BLUE SPECTRA: O II (4639 + 4649 + 4651 + 4662)")
ylabel("BLUE SPECTRA: O II 4642 (plus blended N II, N III)     or    O II (4674 + 4676)")
legend(loc="upper left")
title("O II ratios in V1 multiplet that should be insensitive to (ne, Te)")
gcf().set_size_inches(8, 8)
savefig("oii-insensitive-blue.pdf")

# So it looks like there is definitely some contamination of the 4642 component.  With 4674 + 4676 there is no obvious contamination.  Esteban (2004) have strengths of 0.029, 0.102, 0.015 for N III 4640.64, O II 4641.81, and N II 4643.06, so there is contamination of (0.029 + 0.015)/0.102 = 43%.  With Baldwin (2000), which is a slightly different position, it is (0.0019 + 0.0054)/0.0241 = 30%. 

scatter(totalg, tabg['F4642'], c="r", alpha=0.8)
scatter(totalg, tabg['F4676'], c="b", alpha=0.8)
xlim(0.0, 30.0)
ylim(0.0, 10.0)
xgrid = np.linspace(0.0, 30.0)
plot(xgrid, 0.32*xgrid, "r--", label="O II 4642 theoretical slope: 0.32")
plot(xgrid, 0.32*1.43*xgrid, "r-", label="O II 4642 plus 43% contamination")
plot(xgrid, 0.15*xgrid, "b-", label="O II (4674 + 4676) theoretical slope: 0.15")
xlabel("GREEN SPECTRA: O II (4639 + 4649 + 4651 + 4662)")
ylabel("GREEN SPECTRA: O II 4642 (plus blended N II, N III)     or    O II (4674 + 4676)")
legend(loc="upper left")
title("O II ratios in V1 multiplet that should be insensitive to (ne, Te)")
savefig("oii-insensitive-green.pdf")
gcf().set_size_inches(8, 8)

# So that is the same, but for the green spectra.  It is *far far worse*.  So we are justified in just ignoring the green spectra altogether.

# So once we correct the 4642 component, we now have two density-sensitive O II ratios. 

scatter(1.43*tabb["F4649"]/tabb["F4642"], ratiob, 
        c=radius_b, vmin=30, vmax=100, s=100, alpha=0.6)
xlim(0.7, 1.7)
ylim(0.3, 1.2)
cb = colorbar()
cb.set_label('Radius, arcsec')
title("Both density-sensitive O II V1 ratios")
xlabel("O II 4649 / 4642")
ylabel("O II 4649 / (4639 + 4651 + 4662)")
gcf().set_size_inches(9, 8)
savefig("oii-density-compare.pdf")

# +
from matplotlib.ticker import FormatStrFormatter

plot(radius_b, 1.43*tabb["F4649"]/tabb["F4642"], "ro",
     label="O II 4649 / 4642", alpha=0.5)
plot(radius_b, ratiob, "bo",
     label="O II 4649 / (4639 + 4651 + 4662)", alpha=0.5)
plot(radius_b, 2*ratiob  + 1.43*tabb["F4649"]/tabb["F4642"], "mo",
     label="O II Sum(ratios)", alpha=0.5)
plot(radius_g, clrat, "yo",
     label="[Cl III] 5518 / 5538", alpha=0.3)
xlim(20, 120)
ylim(0.3, 4.0)
title("Density-sensitive ratios versus radius")
xlabel("Radius, arcsec")
ylabel("Ratio")
yscale('log')
xscale('log')
gca().xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
gca().xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
gca().yaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
gca().yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
legend(loc="lower left")
gcf().set_size_inches(8, 12)
savefig("oii-density-vs-radius.pdf")
# -

# The two O II diagnostics both have a very similar fractional spread at each radius.  The [Cl III] diagnostic has a much smaller spread, except for at the smallest radii.  For R > 50'' it looks like it is dominated by the intrinsic spread of densities at each radius.  Since the 4649 / (4639 + 4651 + 4662) ratio has a slightly steeper slope with radius than the 4649 / 4642 ratio, its derived densities will have a slightly smaller spread.  Adding the two ratios (weight of 2:1) doesn;t seem to improve things at all really.

totalg_b = griddata((xg, yg), totalg, (xb, yb), method='linear')
totalb_g = griddata((xb, yb), totalb, (xg, yg), method='linear')
plot(totalb/0.68, totalg_b/0.68, "bo", alpha=0.2, label="Green spectra resampled to Blue positions")
plot(totalb_g/0.68, totalg/0.68, "go", alpha=0.2, label="Blue spectra resampled to Green positions")
plot([0.0, 50.0], [0.0, 50.0], '-r', label="Equal brightness")
xlim(0.0, 40.0)
ylim(0.0, 40.0)
xlabel("BLUE SPECTRA: O II (4639 + 4649 + 4651 + 4662)")
ylabel("GREEN SPECTRA: O II (4639 + 4649 + 4651 + 4662)")
legend(loc="upper left")
title("O II Sum(V1): comparison between Blue and Green spectra")
gcf().set_size_inches(8, 8)

totalg_b = griddata((xg, yg), totalg, (xb, yb), method='linear')
totalb_g = griddata((xb, yb), totalb, (xg, yg), method='linear')
scatter(totalb/0.68, totalg_b/0.68, c=radius_b, vmin=30, vmax=100, s=100, alpha=0.6)
plot([0.0, 50.0], [0.0, 50.0], '-r', label="Equal brightness")
xlim(0.0, 40.0)
ylim(0.0, 40.0)
xlabel("BLUE SPECTRA: O II (4639 + 4649 + 4651 + 4662)")
ylabel("GREEN SPECTRA: O II (4639 + 4649 + 4651 + 4662)")
legend(loc="upper left")
cb = colorbar()
cb.set_label('Radius, arcsec')
title("O II Sum(V1): comparison between Blue and Green spectra")
gcf().set_size_inches(9, 8)

totalg_b = griddata((xg, yg), totalg, (xb, yb), method='linear')
totalb_g = griddata((xb, yb), totalb, (xg, yg), method='linear')
scatter(totalb/0.68, totalb/totalg_b, c=radius_b, vmin=30, vmax=100, s=100, alpha=0.3)
plot([0.0, 50.0], [1.0, 1.0], '-r', label="Equal brightness")
xlim(0.0, 40.0)
ylim(0.0, 2.0)
xlabel("BLUE SPECTRA: O II (4639 + 4649 + 4651 + 4662)")
ylabel("Ratio: BLUE / GREEN")
legend(loc="upper left")
cb = colorbar()
cb.set_label('Radius, arcsec')
title("O II Sum(V1): comparison between Blue and Green spectra")
gcf().set_size_inches(9, 8)
savefig("oii-blue-green-compare.pdf")

# Note the skew towards ratios > 1 for the farthest distances (red points on the left).   This is what is causing the strange "red clump" in the temperature disgnostic plot below. 

hist((totalb/totalg_b)[(radius_b > 60.0) & (radius_b < 85.0)], 
     bins=50, range=[0.0, 2.0], alpha=0.5, normed=True)
hist((totalb/totalg_b)[radius_b < 60.0], 
     bins=50, range=[0.0, 2.0], alpha=0.5, normed=True)
hist((totalb/totalg_b/1.3)[radius_b > 85.0], 
     bins=50, range=[0.0, 2.0], alpha=0.5, normed=True)
xlabel('BLUE / GREEN')
grid()

# This shows that we need to divide the blue spectra by 1.3 for R > 85

totalg_b = griddata((xg, yg), totalg, (xb, yb), method='linear')
totalb_g = griddata((xb, yb), totalb, (xg, yg), method='linear')
mask = radius_b > 85.0
totalb[mask] /= 1.3
scatter(totalb/0.68, totalb/totalg_b, c=radius_b, vmin=30, vmax=100, s=100, alpha=0.3)
plot([0.0, 50.0], [1.0, 1.0], '-r', label="Equal brightness")
xlim(0.0, 40.0)
ylim(0.0, 2.0)
xlabel("BLUE SPECTRA: O II (4639 + 4649 + 4651 + 4662)")
ylabel("Ratio: BLUE / GREEN")
legend(loc="upper left")
cb = colorbar()
cb.set_label('Radius, arcsec')
title("O II Sum(V1): comparison between corrected Blue and Green spectra")
gcf().set_size_inches(9, 8)
savefig("oii-corrected-blue-green-compare.pdf")


# +
def Ratio_v1_4959(T):
    return 6.56e-5*((T/1.e4)**-0.415)*np.exp(29160.0/T)

def Ratio_4363_4959(T):
    return 0.496*np.exp(-32940.0/T)

Ratio_v1_4959(np.array([5.e3, 1.e4, 1.5e4]))

# +
xb, yb, totalb, ratiob = find_xytotrat(tabb)
xg, yg, totalg, ratiog = find_xytotrat(tabg)

b4959 = griddata((xg, yg), tabg["F4959"], (xb, yb), method='linear')
b5007 = griddata((xg, yg), tabg["F5007"], (xb, yb), method='linear')
b4363 = tabb["F4363"]
sumv1bA = totalb + tabb["F4642"]/1.34 + tabb["F4676"]
sumv1bB = totalb/0.68


#plot(b4363/b5007, sumv1bA/b4959, "bo", alpha=0.4)
scatter(b4363/b5007, sumv1bA/b4959, 
        c=radius_b, vmin=30, vmax=100, 
        s=10*sumv1bA, alpha=0.6)
Trange = np.linspace(5000.0, 15000.0, 200)
plot((1./2.918)*Ratio_4363_4959(Trange), Ratio_v1_4959(Trange), 
     "k", label='T(ORL) = T(CEL)')
plot((1./2.918)*Ratio_4363_4959(Trange), Ratio_v1_4959(0.9*Trange), 
     "k--", label='T(ORL) = 0.9 x T(CEL)')
plot((1./2.918)*Ratio_4363_4959(Trange), Ratio_v1_4959(0.8*Trange), 
     "k:", label='T(ORL) = 0.8 x T(CEL)')
xlim(0.002, 0.006)
ylim(0.002, 0.007)
xlabel("[O III] 4363 / 5007")
ylabel("O II Sum(V1) / [O III] 4959")
cb = colorbar()
cb.set_label('Radius, arcsec')
legend(loc="lower right")
title("O II from Blue spectra with [O III] 4959 and 5007 resampled to Blue positions")
gcf().set_size_inches(9, 8)
savefig("oii-oiii-temperature.pdf")
# -

# This is the key figure for the temperature.  It should be compared with Fig 1 and Eq 6 of Peimbert & Peimbert (2013).  *I have now fixed the red points*

#
# This now shows that there is a spread along the curves, which correlates very well with radius.  This is the lock-step variation of the ORL and CEL temperatures. 
#
# There is also a spread across the curves, which is even larger, but that does not correlate with anything that I can see.  This is variation of the ratio of ORL/CEL temperatures, which can be interpreted as t^2 or as the ADF. 
#

# **But we need to de-redden the spectra**   Both are blue/red, so reddening will move us in the same direction as the ADF.   

# +
from matplotlib.ticker import FormatStrFormatter

plot(radius_g, g4363/g5007, "ro",
     label="[O III] 4363/5007", alpha=0.5)
plot(radius_b, sumv1bA/b4959, "bo",
     label="O II Sum(V1) / [O III] 4959", alpha=0.5)
xlim(20, 120)
ylim(0.002, 0.01)
title("Temperature-sensitive ratios versus radius")
xlabel("Radius, arcsec")
ylabel("Ratio")
yscale('log')
xscale('log')
gca().xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
gca().xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
gca().yaxis.set_minor_formatter(FormatStrFormatter("%.3f"))
gca().yaxis.set_major_formatter(FormatStrFormatter("%.3f"))
legend(loc="upper left")
grid()
gcf().set_size_inches(8, 12)
savefig("oii-temperature-vs-radius.pdf")

# +
from matplotlib.ticker import FormatStrFormatter

plot(radius_b, 300*b4363/b5007, "ro",
     label="300 x [O III] 4363/5007", alpha=0.3)
plot(radius_b, b4363/sumv1bA, "bo",
     label="[O III] 4363 / O II Sum(V1)", alpha=0.3)
plot(radius_b, 0.0005*b5007/sumv1bA, "go",
     label="0.0005 x [O III] 5007 / O II Sum(V1)", alpha=0.3)
plot(radius_b, b4649/b4591, "co",
     label="O II V1 4649 / V15 4591", alpha=0.3)
xlim(20, 130)
ylim(0.2, 45.0)
title("Temperature-sensitive ratios versus radius")
xlabel("Radius, arcsec")
ylabel("Ratio")
yscale('log')
xscale('log')
gca().xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
gca().xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
gca().yaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
gca().yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
legend(loc="upper left")
grid()
gcf().set_size_inches(8, 12)
savefig("oii-temperature-three-vs-radius.pdf")
# -

fig, (ax1, ax2) = subplots(1, 2, sharex=True, sharey=True)
scat1 = ax1.scatter(xg, yg, c=g4363/g5007, vmin=0.002, vmax=0.006, s=150, cmap=cm.jet, alpha=0.4)
scat2 = ax2.scatter(xb, yb, c=b4363/b5007, vmin=0.002, vmax=0.006, s=150, cmap=cm.jet, alpha=0.4)
fig.colorbar(scat1, ax=ax1)
fig.colorbar(scat2, ax=ax2)
ax1.set_xlim(0.0, -90.0)
ax1.set_ylim(-90.0, 0.0)
ax1.axis('equal')
ax1.set_title('[O III] 4363/5007 green')
ax2.set_title('[O III] 4363/5007 blue')
ax2.axis('equal')
fig.set_size_inches(12,6)
fig.tight_layout()

fig, (ax1, ax2) = subplots(1, 2, sharex=True, sharey=True)
scat1 = ax1.scatter(xg, yg, c=sumv1bA_g/g4959, vmin=0.002, vmax=0.006, s=150, cmap=cm.jet, alpha=0.4)
scat2 = ax2.scatter(xb, yb, c=sumv1bA/b4959, vmin=0.002, vmax=0.006, s=150, cmap=cm.jet, alpha=0.4)
fig.colorbar(scat1, ax=ax1)
fig.colorbar(scat2, ax=ax2)
ax1.set_xlim(0.0, -90.0)
ax1.set_ylim(-90.0, 0.0)
ax1.axis('equal')
ax1.set_title('O II sum(V1) / [O III] 5007 green')
ax2.set_title('O II sum(V1) / [O III] 5007 blue')
ax2.axis('equal')
fig.set_size_inches(12,6)
fig.tight_layout()

fig, (ax1, ax2) = subplots(1, 2, sharex=True, sharey=True)
scat1 = ax1.scatter(xg, yg, c=g4959, vmin=0.0, s=150, cmap=cm.jet, alpha=0.6)
scat2 = ax2.scatter(xb, yb, c=sumv1bA, vmin=0.0, vmax=30, s=150, cmap=cm.jet, alpha=0.6)
fig.colorbar(scat1, ax=ax1)
fig.colorbar(scat2, ax=ax2)
ax1.set_xlim(0.0, -90.0)
ax1.set_ylim(-90.0, 0.0)
ax1.axis('equal')
ax1.set_title('[O III] 5007 green')
ax2.set_title('O II sum(V1) blue')
ax2.axis('equal')
fig.set_size_inches(12,6)
fig.tight_layout()

fig, (ax1, ax2) = subplots(1, 2, sharex=True, sharey=True)
scat1 = ax1.scatter(xg, yg, c=totalg/0.68, vmin=0.0, vmax=30, s=150, cmap=cm.jet, alpha=0.8)
scat2 = ax2.scatter(xb, yb, c=totalb/0.68, vmin=0.0, vmax=30, s=150, cmap=cm.jet, alpha=0.8)
fig.colorbar(scat1, ax=ax1)
fig.colorbar(scat2, ax=ax2)
ax1.set_xlim(0.0, -90.0)
ax1.set_ylim(-90.0, 0.0)
ax1.axis('equal')
ax1.set_title('O II sum(V1) green')
ax2.set_title('O II sum(V1) blue')
ax2.axis('equal')
fig.set_size_inches(12,6)
fig.tight_layout()

totalg_b = griddata((xg, yg), totalg, (xb, yb), method='nearest')
totalb_g = griddata((xb, yb), totalb, (xg, yg), method='nearest')
plot(totalb/0.68, totalg_b/0.68, "bo", alpha=0.2, label="Green spectra resampled to Blue positions")
plot(totalb_g/0.68, totalg/0.68, "go", alpha=0.2, label="Blue spectra resampled to Green positions")
plot([0.0, 50.0], [0.0, 50.0], '-r', label="Equal brightness")
xlim(0.0, 40.0)
ylim(0.0, 40.0)
xlabel("BLUE SPECTRA: O II (4639 + 4649 + 4651 + 4662)")
ylabel("GREEN SPECTRA: O II (4639 + 4649 + 4651 + 4662)")
legend(loc="upper left")
title("O II Sum(V1): comparison between Blue and Green spectra")
gcf().set_size_inches(8, 8)

# +
g4959 = tabg["F4959"]
g5007 = tabg["F5007"]
b4363 = tabb["F4363"]
sumv1bA = totalb + tabb["F4642"]/1.34 + tabb["F4676"]
sumv1bB = totalb/0.68
g4363 = griddata((xb, yb), tabb["F4363"], (xg, yg), method='linear')
sumv1bA_g = griddata((xb, yb), sumv1bA, (xg, yg), method='linear')
sumv1bB_g = griddata((xb, yb), sumv1bB, (xg, yg), method='linear')

#plot(g4363/g5007, sumv1bA_g/g4959, "bo", alpha=0.4)
scatter(g4363/g5007, sumv1bA_g/g4959,  c=radius_g, vmin=30, vmax=100, s=100, alpha=0.6)
xlim(0.002, 0.008)
ylim(0.002, 0.008)
xlabel("[O III] 4363 / 5007")
ylabel("O II Sum(V1) / [O III] 4959")
title("O II and [O III] 4363 from Blue spectra resampled to Green positions")
gcf().set_size_inches(8, 8)
# -

g4959 = tabg["F4959"]
g5007 = tabg["F5007"]
g4363 = griddata((xb, yb), tabb["F4363"], (xg, yg), method='linear')
sumv1gA = totalg + tabg["F4642"]/1.34 + tabg["F4676"]
sumv1gB = totalg/0.68
#plot(g4363/g5007, sumv1gA/g4959, "bo", alpha=0.4)
scatter(g4363/g5007, sumv1gA/g4959, c=radius_g, vmin=30, vmax=100, s=100, alpha=0.6)
Trange = np.linspace(5000.0, 15000.0, 200)
plot((1./3.)*Ratio_4363_4959(Trange), Ratio_v1_4959(Trange), 
     "k", label='T(ORL) = T(CEL)')
plot((1./3.)*Ratio_4363_4959(Trange), Ratio_v1_4959(0.9*Trange), 
     "k--", label='T(ORL) = 0.9 x T(CEL)')
plot((1./3.)*Ratio_4363_4959(Trange), Ratio_v1_4959(0.8*Trange), 
     "k:", label='T(ORL) = 0.8 x T(CEL)')
xlim(0.002, 0.006)
ylim(0.002, 0.007)
cb = colorbar()
cb.set_label('Radius, arcsec')
legend(loc="lower right")
xlabel("[O III] 4363 / 5007")
ylabel("O II Sum(V1) / [O III] 4959")
title("O II from Green spectra with [O III] 4363 resampled to Green positions")
gcf().set_size_inches(9, 8)

b4959 = griddata((xg, yg), tabg["F4959"], (xb, yb), method='nearest')
b5007 = griddata((xg, yg), tabg["F5007"], (xb, yb), method='nearest')
b4363 = tabb["F4363"]
sumv1bA = totalb + tabb["F4642"]/1.34 + tabb["F4676"]
sumv1bB = totalb/0.68
#plot(b4363/b5007, sumv1bA/b4363, "bo", alpha=0.4)
scatter(b4363/b5007, sumv1bB/b4363, c=radius_b, vmin=30, vmax=100, s=100, alpha=0.6)
xlim(0.002, 0.008)
ylim(0.2, 0.65)
xlabel("[O III] 4363 / 5007")
ylabel("O II Sum(V1) / [O III] 4363")
gcf().set_size_inches(8, 8)

scatter(g4959, g5007)
plot([0.0, 1e4], [0.0, 3.e4])
gcf().set_size_inches(8, 8)

# Now, look at the ratio between multiplets, which should be T-sensitive

b4591 = tabb["F4591"]
scatter(sumv1bA, b4591, c=radius_b, vmin=30, vmax=100, s=50, alpha=0.6) 
xlim(0.0, 35.0)
ylim(0.0, 2.5)
cb = colorbar()
cb.set_label('Radius, arcsec')
xlabel("O II sum(V1)")
ylabel("O II V48a 4591")
gcf().set_size_inches(9, 8)

b4591 = tabb["F4591"]
b4649 = tabb["F4649"]
scatter(b4363/b5007, b4649/b4591, c=radius_b, vmin=30, vmax=100, s=10*b4649, alpha=0.6) 
xlim(0.002, 0.006)
ylim(1.0, 30.0)
yscale('log')
cb = colorbar()
cb.set_label('Radius, arcsec')
xlabel("[O III] 4363 / 5007")
ylabel("O II V1 4649 / V15 4591")
title("BLUE SPECTRA:  O II vs [O III] T diagnostics")
gcf().set_size_inches(9, 8)

b4591 = tabb["F4591"]
b4649 = tabb["F4649"]
scatter(b4363/sumv1bA, b4649/b4591, c=radius_b, vmin=30, vmax=100, s=10*b4649, alpha=0.6) 
xlim(1.5, 5.0)
ylim(0.0, 20.0)
#yscale('log')
cb = colorbar()
cb.set_label('Radius, arcsec')
xlabel("[O III] 4363 / O II sum(V1)")
ylabel("O II V1 4649 / V15 4591")
title("BLUE SPECTRA:  O II vs [O III] T diagnostics")
gcf().set_size_inches(9, 8)

hist(np.log10(g4649/g4591), bins=50, range=[0.0, 1.5], color="y", label="Green spectra")
hist(np.log10(b4649/b4591), bins=50, range=[0.0, 1.5], color="b", label="Blue spectra", alpha=0.6)
xlabel("log 10 (O II V1 4649 / V15 4591)")
legend(loc="upper right")
title("Histogram of O II Temperature diagnostic")

m = b4591 > 1.0

print(*zip(xb[m], yb[m]))

b4591[m][1]

g4591 = tabg["F4591"]
g4649 = tabg["F4649"]
scatter(g4363/g5007, g4649/g4591, c=radius_g, vmin=30, vmax=100, s=10*g4649, alpha=0.6) 
xlim(0.002, 0.006)
ylim(1.0, 30.0)
yscale('log')
cb = colorbar()
cb.set_label('Radius, arcsec')
xlabel("[O III] 4363 / 5007")
ylabel("O II V1 4649 / V48a 4591")
title("GREEN SPECTRA:  O II vs [O III] T diagnostics")
gcf().set_size_inches(9, 8)

hist(sumv1gA/g4959, alpha=0.5, bins=50, range=[0.0, 0.008], color="g")
hist(sumv1bA_g/g4959, alpha=0.5, bins=50, range=[0.0, 0.008], color="b")

hist(sumv1gB, alpha=0.5, bins=50, range=[0.0, 30], color="g")
hist(sumv1bB, alpha=0.5, bins=50, range=[0.0, 30], color="b")

hist(g4959, alpha=0.5, bins=50, range=[0.0, 10000], color="g")
hist(b4959, alpha=0.5, bins=50, range=[0.0, 10000], color="b")

# +
radius_b = np.hypot(xb, yb)
radius_g = np.hypot(xg, yg)

plot(radius_g, g4959, 'ro', alpha=0.4)
plot(radius_b, b4959, 'bo', alpha=0.4)
gcf().set_size_inches(12, 8)
# -

x=np.linspace(2.0, 5.0, 500)
xdata = [2.0, 2.2, 2.4, 2.6, 3.0, 4.0, 4.4, 5.0]
ydata = [0.17, 0.18, 0.21, 0.27, 0.47, 1.17, 1.32, 1.4] 
plot(x, 0.65*(np.tanh(1.2*(x - 3.45)) + 1.2)) # T = 4.0
plot(x, 0.57*(np.tanh(1.2*(x - 3.25)) + 1.36)) # T = 3.6
plot(xdata, ydata, 'o')
grid()

# # Description of data from Manu

# Selected quotes:
#
#
# 1. M42_sr.fits, M42_lr.fits, table_M42_sr.fits y table_M42_lr.fits:  Las dos primeras se corresponden con todos los espectros (unos 8000) observados en el rango rojo, reducidos y calibrados en longitud de onda y flujo (exposición corta y larga, respectivamente). Con estos tienes que tener en cuenta, que los diferentes campos no están reescalados entre sí, ni corregidos de c(Hb).  Esto es así porque yo realizaba eso sobre cada línea de emisión medida.  En cualquier caso, es fácil hacerlo con las tablas asociadas a cada archivo. El mismo posee columnas bastante intuitivas, pero te las explico brevemente:
#
#  * name:  nombre original de la imagen de procedencia
#  * id_ap:  se corresponde con el identificador de la fibra. Lo más importante es que las de cielo van del 401 al 436 y las de campo de la 1 a la 331.
#  * flag: Sirve para no representar fibras solapadas, esto tiene poco sentido en tu caso. Pero creo que pillas la idea. 
#  * dRA, dDEC: Shift en segundos de arco respecto a Ori C
#  * factor1:  Un factor entre exposiciones cortas y largas
#  * factor2: Factor correspondiente a H_alfa en el caso de la exposición corta y a NII6548 en la larga.
#  * c(Hb), err(c(Hb)):  promedio de los valores de la imagen que tú me pasaste de Yusef y O'Dell en la correspondiente posición de fibra.
#
# 2. En el caso de las tablas FITS de cocientes todo es muy parecido. Simplemente que aquí estos cocientes están corregidos de extinción, la nueva columna denominada factor se corresponde con el factor necesario para que todos presenten un mapa suave, pero si lo miras con detalle verás que están cercanos a uno. Y en tu caso al tratarse de una zona del mosaico, probablemente observados en la misma noche, no tendrán cambios.
#
# 3. Y ahora por último, como se trabaja con los datos. Con ambos sería lo mismo (lo primero es lo que no existía y tuve que hacer, para que fuera sencillo obtener todos los espectros de una zona).  La idea es leer las tablas primero, crear una máscara con las condiciones espaciales que uno desea y luego extraer eso de la tabla o de la imagen.  Te adjunto ejemplo.
# 	
#
# He estado pensando varias cosas:
#
# 1. Te recuerdo que tengas cuidado con que los bordes de los espectros no están definidos para todas las fibras. Así que si eliges algún filtro al borde (p. ej. 5755, ten cuidado de que el mismo esté siempre definido)
#
# 2. Rectifico lo que te dije del factor 2.  En principio no importa si todas las fibras que obtienes son de un mismo campo, si no habría que tenerlas en cuenta.  ¿Cómo saber esto?.  Pues en principio lo que debes hacer es usar la propia máscara que obtienes para seleccionar fibras sobre la columna del factor2, si todos los valores son iguales, entonces nada. La importancia de estos factores depende mucho de los cocientes que quieras obtener.
#
#

scatter?
