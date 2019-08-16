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

hdulist = fits.open("Adal-Slits/zorip6rojo_1d.fits")

hdulist.info()

hdu = hdulist[1]

nx, wav0, i0, dwav = [hdu.header[k] for k in ("NAXIS1", "CRVAL1", "CRPIX1", "CD1_1")]
wavs = wav0 + (np.arange(nx) - (i0 - 1))*dwav 

ny, dx = [hdu.header[k] for k in ("NAXIS2", "CD2_2")]
x0, j0 = 0.0, 0.0
xpos = x0 + (np.arange(ny) - (j0 - 1))*dx 


# +
def slit_profile(wav0=5754.6, width=14, trimblue=0, trimred=0, wblue=1.0, wred=1.0, hdu=hdu, wavs=wavs):
    i0 = np.argmin(np.abs(wavs - wav0))
    profile = hdu.data[:, i0-width//2:i0+width//2].mean(axis=1)
    bprofile = hdu.data[:, i0-width+trimblue:i0-width//2].mean(axis=1)
    rprofile = hdu.data[:, i0+width//2:i0+width-trimred].mean(axis=1)
    bgprofile = (wred*rprofile + wblue*bprofile)/(wblue+wred)
    return width*(profile - bgprofile), width*rprofile, width*bprofile

def plot_profile(wav0=5754.6, width=14, trimblue=0, trimred=0, margin=5.0, hdu=hdu, wavs=wavs):
    i0 = np.argmin(np.abs(wavs - wav0))
    wav1, wav2 = wavs[i0-width//2], wavs[i0+width//2]
    wavblue, wavred = wavs[i0-width+trimblue], wavs[i0+width-trimred]
    mask = (wavs >= wavblue - margin) & (wavs <= wavred + margin)
    spec = hdu.data.mean(axis=0)
    plot(wavs[mask], spec[mask])
    ym = spec[mask].max()
    fill_betweenx([1e-20, ym], [wav2, wav2], [wavred, wavred], alpha=0.1, facecolor="r")
    fill_betweenx([1e-20, ym], [wavblue, wavblue], [wav1, wav1], alpha=0.1, facecolor="b")
#    xlim(wavblue, wavred)


# -

plot_profile(5754.6)
grid()

plot_profile(6547.9, trimred=5)
plot_profile(6583.4, width=18, trimblue=3)
plot_profile(6716.4, trimred=4)
plot_profile(6730.8, trimblue=4)
plot(wavs, hdu.data.mean(axis=0), lw=0.2)
yscale('log')
ylim(7.e-16, 8.e-13)
xlim(6450, 6900)
gcf().set_size_inches(12,8)

profile6583, rprofile, bprofile = slit_profile(6583.4, width=18, trimblue=3, wblue=0.0)
plot(xpos, bprofile)
plot(xpos, profile6583, label="6583")
plot(xpos, rprofile)
xlabel('Slit position')
ylabel('Intensity')
grid()
gcf().set_size_inches(12,8)

profile6548, rprofile, bprofile = slit_profile(6547.9, trimred=5, wred=0.0)
plot(xpos, bprofile)
plot(xpos, profile6548, label="6548")
plot(xpos, rprofile)
xlabel('Slit position')
ylabel('Intensity')
grid()
gcf().set_size_inches(12,8)

rat = profile6583/profile6548
plot(xpos, rat)
xlabel('Slit position')
ylabel('Intensity Ratio: 6583/6548')
grid()
ylim(0.0, None)
gcf().set_size_inches(12,8)
print("Mean:", np.mean(rat))
print("Sigma:", np.std(rat), '({:.2f}%)'.format(100*np.std(rat)/np.mean(rat)))

profile5755, rprofile, bprofile = slit_profile(5754.6)
plot(xpos, bprofile)
plot(xpos, profile5755, label="5755")
plot(xpos, rprofile)
xlabel('Slit position')
ylabel('Intensity')
grid()
gcf().set_size_inches(12,8)

plot(xpos, profile5755/profile6583, '.-')
xlabel('Slit position')
ylabel('Intensity Ratio: 5755/6583')
grid()
gcf().set_size_inches(12,8)

profileFeII, rprofile, bprofile = slit_profile(7155.14)
plot(xpos, bprofile)
plot(xpos, 30*profileFeII, label="[Fe II]")
plot(xpos, rprofile)
xlabel('Slit position')
ylabel('Intensity')
grid()
gcf().set_size_inches(12,8)

profile6731, rprofile, bprofile = slit_profile(6730.8, trimblue=4, wblue=0.0)
plot(xpos, bprofile)
plot(xpos, profile6731, label="6731")
plot(xpos, rprofile)
xlabel('Slit position')
ylabel('Intensity')
grid()
gcf().set_size_inches(12,8)

profile6716, rprofile, bprofile = slit_profile(6716.4, trimred=4, wred=0)
plot(xpos, bprofile)
plot(xpos, profile6716, label="6716")
plot(xpos, rprofile)
xlabel('Slit position')
ylabel('Intensity')
grid()
gcf().set_size_inches(12,8)

plot(xpos, profile6716/profile6731)
xlabel('Slit position')
ylabel('Intensity Ratio: 6716/6731')
grid()
gcf().set_size_inches(12,8)

scatter(profile6716/profile6731, profile5755/profile6583)
xlim(0.4, 0.9)
ylim(0.0, 0.05)
xlabel("6716 / 6731")
ylabel("5755 / 6583")
grid()
gcf().set_size_inches(6, 9)

from astropy.table import Table

d = {"Position": xpos, "6583": profile6583, "6548": profile6548, 
     "5755": profile5755, "6716": profile6716, "6731": profile6731}

t = Table(d)

t

t.write("adal-slit6-extract.dat", format="ascii.tab")

t.colnames

ewdata = Table.read("adal-slit6-EW.dat", format="ascii.tab")

plot(xpos, ewdata["W6583"])

Table?

# # Repeat for the blue arm

hdub = fits.open("Adal-Slits/zorip6azul_1d.fits")[1]

nx, wav0, i0, dwav = [hdub.header[k] for k in ("NAXIS1", "CRVAL1", "CRPIX1", "CD1_1")]
bwavs = wav0 + (np.arange(nx) - (i0 - 1))*dwav 

bwavs

plot(bwavs, hdub.data.mean(axis=0), "b", label="zorip6azul")
plot(wavs, hdu.data.mean(axis=0), "r", label="zorip6rojo")
yscale("log")
#xscale("log")
xlim(4200.0, 8600.0)
ylim(1e-16, 1e-12)
xticks(np.arange(4500, 9000, 500))
minorticks_on()
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
grid()
legend()
gcf().set_size_inches(15, 6)
tight_layout()
savefig("Adal_Slit6_full_spectrum.pdf")

plot(wavs, hdu.data.mean(axis=0), "r", label="zorip6rojo")
#yscale("log")
#xscale("log")
xlim(6600.0, 8000.0)
ylim(0, 0.4e-14)
xticks(np.arange(6600, 8000, 100))
minorticks_on()
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
grid()
legend()
gcf().set_size_inches(15, 6)
tight_layout()

plot(bwavs, hdub.data.mean(axis=0), "b", label="zorip6azul")
xlim(4150, 4850)
ylim(0.0, 2e-15)
minorticks_on()
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
grid(ls='-', c='b', lw=0.6, alpha=0.3)
grid(ls='-', c='b', lw=0.6, alpha=0.05, which='minor')
gcf().set_size_inches(12, 8)
tight_layout()

plot(bwavs, hdub.data.mean(axis=0), "b", label="zorip6azul")
xlim(4580, 4740)
ylim(0.3e-15, 0.8e-15)
minorticks_on()
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
grid(ls='-', c='b', lw=0.6, alpha=0.3)
grid(ls='-', c='b', lw=0.6, alpha=0.05, which='minor')
gcf().set_size_inches(12, 8)
tight_layout()

hdub.data.shape

mask = (bwavs > 4652) & (bwavs < 4655)
nblock = 10
spread = 1.0
ny, nx = hdub.data.shape
nspec = ny//nblock
for i in range(nspec):
    spec = hdub.data[nblock*i:nblock*(i+1),:].mean(axis=0)
    spec /= spec[mask].mean()
    if spec[mask].std() > 0.3:
        continue
#    print(i, spec[mask].std())
    plot(bwavs, spec+spread*i)
xlim(4880, 4980)
ylim(0.0, 5.0 + spread*nspec)
minorticks_on()
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
grid(ls='-', c='b', lw=0.6, alpha=0.3)
grid(ls='-', c='b', lw=0.6, alpha=0.05, which='minor')
gcf().set_size_inches(12, 8)
tight_layout()
gcf().savefig("Adal-slit6-multispec-4880-4980.pdf")

mask = (bwavs > 4652) & (bwavs < 4655)
nblock = 10
spread = 1.0
ny, nx = hdub.data.shape
nspec = ny//nblock
for i in range(nspec):
    spec = hdub.data[nblock*i:nblock*(i+1),:].mean(axis=0)
    spec /= spec[mask].mean()
    if spec[mask].std() > 0.3:
        continue
#    print(i, spec[mask].std())
    plot(bwavs, spec+spread*i)
xlim(4780, 4880)
ylim(0.0, 5.0 + spread*nspec)
minorticks_on()
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
grid(ls='-', c='b', lw=0.6, alpha=0.3)
grid(ls='-', c='b', lw=0.6, alpha=0.05, which='minor')
gcf().set_size_inches(12, 8)
tight_layout()
gcf().savefig("Adal-slit6-multispec-4780-4880.pdf")

mask = (bwavs > 4652) & (bwavs < 4655)
nblock = 10
spread = 1.0
ny, nx = hdub.data.shape
nspec = ny//nblock
for i in range(nspec):
    spec = hdub.data[nblock*i:nblock*(i+1),:].mean(axis=0)
    spec /= spec[mask].mean()
    if spec[mask].std() > 0.3:
        continue
#    print(i, spec[mask].std())
    plot(bwavs, spec+spread*i)
xlim(4680, 4780)
ylim(0.0, 5.0 + spread*nspec)
minorticks_on()
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
grid(ls='-', c='b', lw=0.6, alpha=0.3)
grid(ls='-', c='b', lw=0.6, alpha=0.05, which='minor')
gcf().set_size_inches(12, 8)
tight_layout()
gcf().savefig("Adal-slit6-multispec-4680-4780.pdf")

mask = (bwavs > 4652) & (bwavs < 4655)
nblock = 10
spread = 1.0
ny, nx = hdub.data.shape
nspec = ny//nblock
for i in range(nspec):
    spec = hdub.data[nblock*i:nblock*(i+1),:].mean(axis=0)
    spec /= spec[mask].mean()
    if spec[mask].std() > 0.3:
        continue
#    print(i, spec[mask].std())
    plot(bwavs, spec+spread*i)
xlim(4580, 4680)
ylim(0.0, 5.0 + spread*nspec)
minorticks_on()
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
grid(ls='-', c='b', lw=0.6, alpha=0.3)
grid(ls='-', c='b', lw=0.6, alpha=0.05, which='minor')
gcf().set_size_inches(12, 8)
tight_layout()
gcf().savefig("Adal-slit6-multispec-4580-4680.pdf")

mask = (bwavs > 4652) & (bwavs < 4655)
nblock = 10
spread = 0.4
ny, nx = hdub.data.shape
nspec = ny//nblock
for i in range(nspec):
    spec = hdub.data[nblock*i:nblock*(i+1),:].mean(axis=0)
    spec /= spec[mask].mean()
    if spec[mask].std() > 0.3:
        continue
#    print(i, spec[mask].std())
    plot(bwavs, spec+spread*i)
xlim(4480, 4580)
ylim(0.0, 5.0 + spread*nspec)
minorticks_on()
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
grid(ls='-', c='b', lw=0.6, alpha=0.3)
grid(ls='-', c='b', lw=0.6, alpha=0.05, which='minor')
gcf().set_size_inches(12, 8)
tight_layout()
gcf().savefig("Adal-slit6-multispec-4480-4580.pdf")

# This shows the He II absorption at 4541 Ang.

mask = (bwavs > 4652) & (bwavs < 4655)
nblock = 10
spread = 1.0
ny, nx = hdub.data.shape
nspec = ny//nblock
for i in range(nspec):
    spec = hdub.data[nblock*i:nblock*(i+1),:].mean(axis=0)
    spec /= spec[mask].mean()
    if spec[mask].std() > 0.3:
        continue
#    print(i, spec[mask].std())
    plot(bwavs, spec+spread*i)
xlim(4380, 4480)
ylim(0.0, 5.0 + spread*nspec)
minorticks_on()
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
grid(ls='-', c='b', lw=0.6, alpha=0.3)
grid(ls='-', c='b', lw=0.6, alpha=0.05, which='minor')
gcf().set_size_inches(12, 8)
tight_layout()
gcf().savefig("Adal-slit6-multispec-4380-4480.pdf")

mask = (bwavs > 4652) & (bwavs < 4655)
nblock = 10
spread = 10.0
ny, nx = hdub.data.shape
nspec = ny//nblock
for i in range(nspec):
    spec = hdub.data[nblock*i:nblock*(i+1),:].mean(axis=0)
    spec /= spec[mask].mean()
    if spec[mask].std() > 0.3:
        continue
#    print(i, spec[mask].std())
    plot(bwavs, spec+spread*i)
xlim(4280, 4380)
ylim(0.0, 5.0 + spread*nspec)
minorticks_on()
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
grid(ls='-', c='b', lw=0.6, alpha=0.3)
grid(ls='-', c='b', lw=0.6, alpha=0.05, which='minor')
gcf().set_size_inches(12, 8)
tight_layout()
gcf().savefig("Adal-slit6-multispec-4280-4380.pdf")

mask = (bwavs > 4652) & (bwavs < 4655)
nblock = 10
spread = 1.0
ny, nx = hdub.data.shape
nspec = ny//nblock
for i in range(nspec):
    spec = hdub.data[nblock*i:nblock*(i+1),:].mean(axis=0)
    spec /= spec[mask].mean()
    if spec[mask].std() > 0.3:
        continue
#    print(i, spec[mask].std())
    plot(bwavs, spec+spread*i)
xlim(4180, 4280)
ylim(0.0, 5.0 + spread*nspec)
minorticks_on()
xlabel("Wavelength, Angstrom")
ylabel("Flux, erg/cm2/s/A")
grid(ls='-', c='b', lw=0.6, alpha=0.3)
grid(ls='-', c='b', lw=0.6, alpha=0.05, which='minor')
gcf().set_size_inches(12, 8)
tight_layout()
gcf().savefig("Adal-slit6-multispec-4180-4280.pdf")

hdu = fits.open("Adal-Slits/zorip6azul_1d.fits")[1]
#hdu.data = hdu.data[:120,:]
nx, wav0, i0, dwav = [hdub.header[k] for k in ("NAXIS1", "CRVAL1", "CRPIX1", "CD1_1")]
wavs = wav0 + (np.arange(nx) - (i0 - 1))*dwav 
ny, dx = [hdu.header[k] for k in ("NAXIS2", "CD2_2")]
x0, j0 = 0.0, 0.0
xpos = x0 + (np.arange(ny) - (j0 - 1))*dx 
#xpos = xpos[:120]

hdur = fits.open("Adal-Slits/zorip6rojo_1d.fits")[1]
#hdu.data = hdu.data[:120,:]
nx, wav0, i0, dwav = [hdur.header[k] for k in ("NAXIS1", "CRVAL1", "CRPIX1", "CD1_1")]
rwavs = wav0 + (np.arange(nx) - (i0 - 1))*dwav 
ny, dx = [hdu.header[k] for k in ("NAXIS2", "CD2_2")]
x0, j0 = 0.0, 0.0
xpos = x0 + (np.arange(ny) - (j0 - 1))*dx 
#xpos = xpos[:120]

hdu5 = fits.open("Adal-Slits/zorip5azul_1d.fits")[1]
ny, dx = [hdu.header[k] for k in ("NAXIS2", "CD2_2")]
x0, j0 = 0.0, 0.0
xpos5 = x0 + (np.arange(ny) - (j0 - 1))*dx 


plot_profile(4590.97, width=10)
plot_profile(4590.97, width=10, hdu=hdu5)
title("O II V15 4591 - Temperature indicator")

plot_profile(4649.13, width=8, trimred=3)
plot_profile(4649.13, width=8, trimred=3, hdu=hdu5)
title("O II V1 4649 - Density indicator")

plot_profile(4650.84, width=8, trimblue=3)
title("O II V1 4651")

plot_profile(4661.63, width=8)
title("O II V1 4662")

plot_profile(4641.6, width=12, trimblue=2)
title("O II V1 4642 plus N III on the flanks")

plot_profile(4638.86, width=8, trimred=2)
title("O II V1 4639")

plot_profile(4676.24, width=8)
title("O II V1 4676")

plot_profile(4673.7, width=8, trimred=2)
title("O II V1 4674")

plot_profile(4696, width=8, margin=2)
title("O II V1 4696")

# +
profile4639, _, _ = slit_profile(4638.86, width=8, trimred=2)
profile4642, _, _ = slit_profile(4641.6, width=12, trimblue=2)
profile4649, _, _ = slit_profile(4649.13, width=9, trimred=4)
profile4651, _, _ = slit_profile(4650.84, width=8, trimblue=3)
profile4662, _, _ = slit_profile(4661.63, width=8)
profile4676, _, _ = slit_profile(4676.24, width=8)
profile4674, _, _ = slit_profile(4673.7, width=8, trimred=2)

denratio = profile4649 / (profile4639 + profile4651 + profile4662)
denratio2 = profile4649 / profile4642
sum_4650 = profile4639 + profile4642 + profile4649 + profile4651 + profile4662 + profile4676 + profile4674

# -

plot(denratio2, denratio, 'o')
xlim(0.0, 2.5)
ylim(0.0, 2.5)
title("Both density-sensitive O II V1 ratios")
xlabel("O II 4649 / 4642")
ylabel("O II 4649 / (4639 + 4651 + 4662)")
gcf().set_size_inches(6, 6)

xgrid = np.linspace(0.0, 1.6e-14)
insensitive_sum = profile4639 + profile4649 + profile4651 + profile4662
scatter(insensitive_sum, profile4642, c='r', alpha=0.8)
plot(xgrid, 0.32*xgrid, "r", label="O II 4642 theoretical slope: 0.32")
scatter(insensitive_sum, profile4676, c='b', alpha=0.8)
plot(xgrid, 0.13*xgrid, "b", label="O II 4676 theoretical slope: 0.13")
scatter(insensitive_sum, profile4674, c='y', alpha=0.8)
plot(xgrid, 0.02*xgrid, "y", label="O II 4674 theoretical slope: 0.02")
xlim(0.0, 1.6e-14)
ylim(0.0, 0.6e-14)
xlabel("O II (4639 + 4649 + 4651 + 4662)")
ylabel("O II 4642  or  4674  or  4676")
legend(loc="upper left")
title("O II ratios in V1 multiplet that should be insensitive to (ne, Te)")
gcf().set_size_inches(8, 8)
savefig("oii-insensitive-adal-slit6.pdf")

# There is evidence that the 4642 intensity is higher than the theoretical prediction.  But that line still has problems with the adjacent N III.  I need to fit Gaussians to do this properly. 

# The 4676 line seems slightly below the theoretical slope for the higher intensities.  Not sure if this is significant. 

# +
plot(xpos, profile4639, label="4639")
plot(xpos, profile4642, label="4642")
plot(xpos, profile4651, label="4651")
plot(xpos, profile4662, label="4662")
plot(xpos, profile4676, label="4676")
plot(xpos, profile4674, label="4674")
plot(xpos, profile4649, "k", label="4649")
plot(xpos, sum_4650, "k", lw=2, label="sum")

xlabel('Slit position')
ylabel('Intensity')
ylim(-1e-15, 2.0e-14)
legend()
grid()
gcf().set_size_inches(12,8)
# -

plot_profile(4740, width=8)
title("[Ar IV] 4740")

plot_profile(4711, width=8, trimred=3)
title("[Ar IV] 4711")

# +
profile4711, _, _ = slit_profile(4711, width=8, trimred=3)
profile4740, _, _ = slit_profile(4740, width=8)
ariv_ratio = profile4711/profile4740

scatter(ariv_ratio, denratio, c=xpos, s=100, alpha=0.6)
colorbar()
xlim(0.0, 2.0)
ylim(0.0, 2.0)
# -

nb = 14
avArIV = block_average(ariv_ratio, nb)
avOII = block_average(denratio, nb)
avx = block_average(xpos, nb)
scatter(avOII, avArIV, c=avx, s=100, alpha=0.6)
plot(avOII, avArIV, '-k')
colorbar()
xlabel("O II 4649 / (4639 + 4651 + 4662)")
ylabel("[Ar IV] 4711/4740")
xlim(0.0, 2.0)
ylim(0.0, 2.0)

profile4642, rprofile, bprofile = slit_profile(4711, width=8, trimred=3, hdu=hdu5)
plot(xpos, bprofile)
plot(xpos, profile4642, label="4711")
plot(xpos, rprofile)
plot(xpos, np.abs(bprofile - rprofile))
xlabel('Slit position')
ylabel('Intensity')
ylim(-1e-14, 2e-14)
grid()
gcf().set_size_inches(12,8)

profile4642, rprofile, bprofile = slit_profile(4740, width=8, hdu=hdu5)
plot(xpos, bprofile)
plot(xpos, profile4642, label="4740")
plot(xpos, rprofile)
plot(xpos, np.abs(bprofile - rprofile))
xlabel('Slit position')
ylabel('Intensity')
ylim(-1e-14, 2e-14)
grid()
gcf().set_size_inches(12,8)

plot_profile(5876, width=14, hdu=hdur, wavs=rwavs)
title("He I 5876")

plot_profile(6678, width=14, hdu=hdur, wavs=rwavs)
title("He I 6678")

profile6678, rprofile, bprofile = slit_profile(6678, width=14, hdu=hdur, wavs=rwavs)
plot(xpos, bprofile)
plot(xpos, profile6678, label="6678")
plot(xpos, rprofile)
plot(xpos, np.abs(bprofile - rprofile))
xlabel('Slit position')
ylabel('Intensity')
ylim(-1e-14, 2e-13)
grid()
gcf().set_size_inches(12,8)

profile5876, rprofile, bprofile = slit_profile(5876, width=14, hdu=hdur, wavs=rwavs)
plot(xpos, bprofile)
plot(xpos, profile5876, label="5876")
plot(xpos, rprofile)
plot(xpos, np.abs(bprofile - rprofile))
xlabel('Slit position')
ylabel('Intensity')
ylim(-1e-14, 5e-13)
grid()
gcf().set_size_inches(12,8)

hei_ratio = profile6678 / profile5876
plot(xpos, hei_ratio)

scatter(hei_ratio, denratio, c=xpos, s=100, alpha=0.6)
colorbar()
xlim(0.28, 0.36)
ylim(0.0, 2.0)

avHeI = block_average(hei_ratio, nb)
scatter(avHeI, avOII, c=avx, s=100, alpha=0.6)
plot(avHeI, avOII, '-k')
colorbar()
xlim(0.28, 0.36)
ylim(0.0, 2.0)
xlabel('He I 6678 / 5876')
ylabel("O II 4649 / (4639 + 4651 + 4662)")
title("O II density vs He I density")

# The He I ratio will go down when de-reddened.  If chb ~ 1 then by 10^0.1 = 1.25 => 0.25.  This is exactly the expected ratio for 1e4 pcc, increasing by 10% or so as density drops to 1e3 pcc. 

profile4642, rprofile, bprofile = slit_profile(4641.6, width=12, trimblue=2)
plot(xpos, bprofile)
plot(xpos, profile4642, label="4642")
plot(xpos, rprofile)
plot(xpos, np.abs(bprofile - rprofile))
xlabel('Slit position')
ylabel('Intensity')
ylim(-1e-14, 2e-14)
grid()
gcf().set_size_inches(12,8)

profile4591, rprofile, bprofile = slit_profile(4590.97, width=10)
plot(xpos, bprofile)
plot(xpos, profile4591, label="4591")
plot(xpos, rprofile)
plot(xpos, np.abs(bprofile - rprofile))
xlabel('Slit position')
ylabel('Intensity')
ylim(-1e-14, 2e-14)
grid()
gcf().set_size_inches(12,8)

profile4649, rprofile, bprofile = slit_profile(4649.13, width=9, trimred=4)
plot(xpos, bprofile)
plot(xpos, profile4649, label="4649")
plot(xpos, rprofile)
plot(xpos, np.abs(bprofile - rprofile))
xlabel('Slit position')
ylabel('Intensity')
ylim(-1e-14, 2e-14)
grid()
gcf().set_size_inches(12,8)


def block_average(x, n):
    m = len(x)
    return x.reshape((m/n, n)).mean(axis=-1)


print(block_average(xpos, 7))

#mask = rprofile > 2e-14
#profile4649[mask] = profile4649[~mask].mean()
#profile4591[mask] = profile4591[~mask].mean()
x = block_average(xpos, 14)
r = block_average(profile4649, 14)/block_average(profile4591, 14)
plot(xpos, profile4649/profile4591, alpha=0.1)
plot(x, r)
ylim(0.0, 15.0)

plot_profile(5006.6, width=12)

plot_profile(4958.5, width=12)

plot_profile(4363.21, width=10)

plot_profile(4931.32 - 0.5, width=12)

# This is the weakest [O III] nebular line at 4931.32 (or maybe 4931.227), but unfortunately it is blended with an [Fe III] line at 4930.54, so we can't really use it. 

# In principle the ratios 4931:4959:5007 should be:

A4931 = 2.41e-06 
A4959 = 6.21e-03 + 4.57e-06
A5007 = 1.81e-02 + 3.52e-05
print("4931:4959:5007 = {:.5f}:1.00000:{:.5f}".format(A4931/A4959, A5007/A4959))

profile4931, rprofile, bprofile = slit_profile(4931.32 - 0.5, width=12)
plot(xpos, bprofile)
plot(xpos, profile4931, label="4931")
plot(xpos, rprofile)
plot(xpos, np.abs(bprofile - rprofile))
xlabel('Slit position')
ylabel('Intensity')
ylim(-1e-15, 2e-14)
grid()
legend()
gcf().set_size_inches(12,8)

plot(xpos, profile4931/profile4959)
axhline(0.0)
axhline(0.00039)
xlabel("Position")
ylabel("([O III] 4931 + [Fe III] 4931) / [O III] 4959") 

# The horizontal line is the theoretical [O III] 4931/4959 ratio.  The excess is presumably due to the [Fe III] line. 

profile4363, rprofile, bprofile = slit_profile(4363.21, width=10)
plot(xpos, bprofile)
plot(xpos, profile4363, label="4363")
plot(xpos, rprofile)
plot(xpos, np.abs(bprofile - rprofile))
xlabel('Slit position')
ylabel('Intensity')
#ylim(-1e-14, 2e-14)
grid()
legend()
gcf().set_size_inches(12,8)

plot(xpos, denratio)
plot(xpos, denratio2)
ylim(0.0, None)

plot(xpos, profile4363/sum_4650)
ylim(0.0, None)

plot(xpos, 300*profile4649, label="O II 4649 x 300")
plot(xpos, 0.6*profile4959, label="[O III] 4959 x 0.6")
#plot(xpos, 180*sum_4650, label="sum(v1)")
legend()
xlabel("Position")
ylabel("[O III] and O II profiles")
ylim(0.0, None)

# This is a recreation of Adal's Fig 5(e) from Mesa Delgado et al (2008). 

# It looks very similar, apart from the normalization. 

plot(xpos, 0.001*profile4959/profile4649)
xlabel("Position")
ylabel("[O III] 4959 / 1000 x O II 4649")
ylim(0.55, 1.4)
gcf().set_size_inches(12,3)

# So this is a recreation of the bottom panel of Adal's Fig 5(e).  It mainly resembles it, but mine is much noisier.  The average value of 0.7 to 0.9 is in agreement. This implies that 4649/4959 is about 1.25e-3. 

plot(xpos, profile4649/profile4959, label="V1 4649")
plot(xpos, sum_4650/profile4959, label="Sum(V1)")
xlabel("Position")
legend()
ylabel("O II lines / [O III] 4959")
axhline(0.003440, ls='--', c='g')
axhline(0.0011112, ls='--', c='b')
ylim(0.0, 0.006)

# So the above plot shows this, together with the sum of all the V1 components.  The dashed lines show the corresponding results from Esteban.  Yes - it now works!!
#

plot(xpos, profile4649/profile4959/0.0011112, label="V1 4649 / average")
plot(xpos, sum_4650/profile4959/0.003440, label="Sum(V1) / average")
xlabel("Position")
legend()
ylabel("O II lines / [O III] 4959")
ylim(0.0, 2.0)

# The above is the same as the previous graph, but normalizing both ratios to the Esteban values.  The Sum(V1) is slightly less noisy in most parts but has regions where it has peaks.
#

profile4959, rprofile, bprofile = slit_profile(4958.5, width=12)
profile5007, rprofile, bprofile = slit_profile(5006.6, width=12)
#plot(xpos, profile4363/profile5007)
plot(xpos, profile4363/(2.91817*profile4959))
plot(xpos, 5e8*profile4959)
ylim(0.00, 0.0045)

scatter(profile4959, profile4363/(2.91817*profile4959), c=xpos, 
        s=100, alpha=0.6)
xlim(0.0, 7e-12)
xlabel('[O III] 4959')
ylabel('[O III] 4363 / 5007')
colorbar()
ylim(0.0025, 0.0045)

# The faint parts of the nebula have the lowest 4363/5007 ratio.  The bright parts bifurcate into two regions.  The region directly S of th1c has a slightly elevated ratio.  The region around Orion S has a more highly elevated ratio. 

scatter(profile4959, sum_4650/profile4959, c=xpos, 
        s=100, alpha=0.6)
xlim(0.0, 7e-12)
xlabel('[O III] 4959')
ylabel('O II Sum(V1) / [O III] 4959')
colorbar()
ylim(0.0025, 0.0045)

# Now we repeat, but using the O II / [O III] T-sensitive ratio.  This shows similar behaviour to the previous plot for the bright parts of the slit: bifurcation into two zones, except upside down because the ratio falls with T.  For the faint parts of the slit the data are too noisy to say anything.

scatter(profile4959, profile4649/profile4591, c=xpos, 
        s=100, alpha=0.6)
xlim(0.0, 7e-12)
xlabel('[O III] 4959')
ylabel('O II V1 4649 / O II V15 5491')
colorbar()
ylim(0, 15)

# And again with the pure O II T-sensitive ratio.   The same effect is seen, but a bit noisier.  Note, however, that according to Fig. 21 of Fang (2013), the absolute T varies between 8000 K (ratio = 3) 4500 K (ratio = 6) and 3000 K (ratio = 9). 

scatter(profile4959, profile4591/profile4959, c=xpos, 
        s=100, alpha=0.6)
xlim(0.0, 7e-12)
colorbar()
ylim(0, 0.0005)

scatter(profile4363/(2.91817*profile4959), 
        profile4591/profile4959, c=xpos, 
        s=100, alpha=0.6)
xlim(0.0025, 0.0045)
colorbar()
ylim(0, 0.0005)

# So the 4591 O II recombination line does not show any significant variations with respect to the [O III] 4959 line.  This implies that they have a similar T dependence.  Either that, or there is a conicidental cancellation. 

# ## Repeat the above, but binning the positions

nb = 14
av4959 = block_average(profile4959, nb)
avV1 = block_average(sum_4650, nb)
avx = block_average(xpos, nb)
scatter(av4959, avV1/av4959, c=avx, s=100, alpha=0.6)
plot(av4959, avV1/av4959, '-k')
xlim(0.0, 7e-12)
xlabel('[O III] 4959')
ylabel('O II Sum(V1) / [O III] 4959')
colorbar()
ylim(0.00, 0.0045)
title('O++: R(ORL–CEL) vs S(CEL)')

av4649 = block_average(profile4649, nb)
av4591 = block_average(profile4591, nb)
scatter(av4959, av4649/av4591, c=avx, 
        s=100, alpha=0.6)
plot(av4959, av4649/av4591, '-k')
xlim(0.0, 7e-12)
xlabel('[O III] 4959')
ylabel('O II V1 4649 / O II V15 4591')
colorbar()
ylim(0, 15)
title('O++: R(ORL) vs S(CEL)')

av4363 = block_average(profile4363, nb)
scatter(av4959, av4363/(2.91817*av4959), c=avx, 
        s=100, alpha=0.6)
plot(av4959, av4363/(2.91817*av4959), '-k')
xlim(0.0, 7e-12)
xlabel('[O III] 4959')
ylabel('[O III] 4363 / 5007')
colorbar()
ylim(0.0025, 0.0045)
title('O++: R(CEL) vs S(CEL)')

scatter(av4363/(2.91817*av4959), av4649/av4591, c=avx, 
        s=100, alpha=0.6)
plot(av4363/(2.91817*av4959), av4649/av4591, '-k')
xlim(0.0025, 0.0045)
ylabel('O II V1 4649 / O II V15 4591')
xlabel('[O III] 4363 / 5007')
colorbar()
ylim(0, 15)
title('O++: R(ORL) vs R(CEL)')


# +
def Ratio_v1_4959(T):
    return 6.56e-5*((T/1.e4)**-0.415)*np.exp(29160.0/T)

def Ratio_4363_4959(T):
    return 0.496*np.exp(-32940.0/T)

Ratio_v1_4959(np.array([5.e3, 1.e4, 1.5e4]))
# -

scatter(av4363/(2.91817*av4959), avV1/av4959, c=avx, 
        s=100, alpha=0.6)
plot(av4363/(2.91817*av4959), avV1/av4959, '-k')
xlim(0.0025, 0.0045)
ylabel('O II Sum(V1) / [O III] 4959')
xlabel('[O III] 4363 / 5007')
colorbar()
ylim(0.0025, 0.0045)
Trange = np.linspace(5000.0, 15000.0, 200)
plot((1./2.918)*Ratio_4363_4959(Trange), Ratio_v1_4959(Trange), 
     "k", label='T(ORL) = T(CEL)')
plot((1./2.918)*Ratio_4363_4959(Trange), Ratio_v1_4959(0.95*Trange), 
     "k:", label='T(ORL) = 0.95 x T(CEL)')
plot((1./2.918)*Ratio_4363_4959(Trange), Ratio_v1_4959(0.9*Trange), 
     "k--", label='T(ORL) = 0.9 x T(CEL)')
plot((1./2.918)*Ratio_4363_4959(Trange), Ratio_v1_4959(0.85*Trange), 
     "k-.", label='T(ORL) = 0.85 x T(CEL)')
axvline((1./2.918)*Ratio_4363_4959(9000), alpha=0.4)
axvline((1./2.918)*Ratio_4363_4959(8500), alpha=0.4)
axvline((1./2.918)*Ratio_4363_4959(8000), alpha=0.4)
title('O++: R(ORL–CEL) vs R(CEL)')

# So that last one may need to be modified by reddening.



# ### Older stuff

plot(xpos, profile5007/profile4959)
ylim(0.0, 3.1)

plot(xpos, profile4959/sum_4650)
ylim(0.0, None)

scatter(profile4363/sum_4650, profile4649/profile4591, 
        s=5e15*sum_4650, alpha=0.5)
ylim(0.0, 15)

scatter(profile4363/(2.918*profile4959), profile4649/profile4591, 
        s=1e15*profile4363, alpha=0.5)
ylim(0.0, 15)
xlim(0.0025, 0.0045)

# +
scatter(profile4363/(2.918*profile4959), sum_4650/profile4959, 
        c=xpos, #cmap=cm.hsv,
        s=2e15*profile4363, alpha=0.6)
Trange = np.linspace(5000.0, 15000.0, 200)
plot((1./2.918)*Ratio_4363_4959(Trange), Ratio_v1_4959(Trange), 
     "k", label='T(ORL) = T(CEL)')
plot((1./2.918)*Ratio_4363_4959(Trange), Ratio_v1_4959(0.9*Trange), 
     "k--", label='T(ORL) = 0.9 x T(CEL)')
plot((1./2.918)*Ratio_4363_4959(Trange), Ratio_v1_4959(0.8*Trange), 
     "k:", label='T(ORL) = 0.8 x T(CEL)')
axvline((1./2.918)*Ratio_4363_4959(9000))
axvline((1./2.918)*Ratio_4363_4959(8000))
xlabel("[O III] 4363 / 5007")
ylabel("O II Sum(V1) / [O III] 4959")
colorbar()
ylim(0.002, 0.007)
xlim(0.002, 0.006)

legend(loc="center right")
title("O II and [O III] from Adal's spectra")
gcf().set_size_inches(9, 8)
savefig("oii-oiii-temperature-adal-slit6.pdf")
# -

split?

# So this looks fine now.

# Check against the values in Esteban (before de-reddening):

est_sumv1 = 0.053 + 0.096 + 0.146 + 0.049 + 0.064 + 0.011 + 0.033
est_4959 = 131.389
est_5007 = 398.147
est_4363 = 1.129
print("Esteban 4363 / 5007 =", est_4363/est_5007)
print("Esteban V1 / 4959 =", est_sumv1/est_4959)
print("Esteban 4649 / 4959 =", 0.146/est_4959)

# Partition the points into five sets based on position.  This is motivated by the graphs above that seem to show 4 regimes in the T-sensitive ratios, with the fifth being where there is a star.

splits = [40, 50, 62, 67, 80, 95]
specs = [section.mean(axis=0) for section in split(hdu.data, splits)]
del(specs[3])
xx = [section.mean()/xpos.max() for section in split(xpos, splits)]
del(xx[3])

xx

# +
contmask = (wavs > 4545) & (wavs < 4560)
normmask = (wavs > 4955) & (wavs < 4962)

def mycmap(x):
    xmin, xmax = 0.15, 1.0
    return cm.RdPu(xmin + x*(xmax-xmin))

nplots = len(specs)
fig, axes = subplots(nplots, 1)
for i, spec in enumerate(specs):
    for ax in axes:
        ax.plot(wavs + 0.4*(1+i-nplots/2), 
                (spec - spec[contmask].mean())/spec[normmask].max(),
                c=mycmap((xx[i]-min(xx))/(max(xx) - min(xx))), 
#                c=cm.jet(xx[i]), 
                lw=4-0.5*i, 
                alpha=0.7)

scale = 2e-2
axes[0].set_ylim(-0.1*scale, scale)
axes[0].set_xlim(4350, 4380)
scale = 3e-4
axes[1].set_ylim(-0.1*scale, scale)
axes[1].set_xlim(4587, 4617)
scale = 1.5e-3
axes[2].set_ylim(-0.1*scale, scale)
axes[2].set_xlim(4617, 4647)
scale = 1.5e-3
axes[3].set_ylim(-0.1*scale, scale)
axes[3].set_xlim(4645, 4680)
scale = 1.5e-3
axes[4].set_ylim(-0.1*scale, scale)
axes[4].set_xlim(4708, 4743)
scale = 1.1
axes[5].set_ylim(-0.1*scale, scale)
axes[5].set_xlim(4950, 4970)

gcf().set_size_inches(12, 20)
tight_layout()
savefig("adal-slit6-oii-multispec.pdf")
# -

array_split?

# # Profiles through the WFC3 filter passbands

from wfc3_utils import get_filter, get_interpolated_filter


import numpy as np

# Get all the filters that fall in the red spectrum, interpolated onto Adal's wavelengths.

wfc3_filters = ["F673N", "FQ672N", "FQ674N", "FQ575N", "F547M",
                "F658N", "F656N"]
T = {f: get_interpolated_filter(f, wavs) for f in wfc3_filters}


# Plot them together with the integrated spectrum from the entire slit.

# +
def Tplot(fname, wavs=wavs):
    plot(wavs, T[fname], lw=3, alpha=0.6, label=fname)

Tplot("F673N")
Tplot("FQ672N")
Tplot("FQ674N")
fill_between(wavs, 2.e13*hdu.data.mean(axis=0), np.zeros_like(wavs),
             alpha=0.4, facecolor="k")
xlim(6600, 6900)
ylim(0.0, 0.3)
xlabel("Wavelength, Angstrom")
ylabel("Transmission")
legend()
gcf().set_size_inches(10, 8)
# -

Tplot("F656N")
Tplot("F658N")
fill_between(wavs, 0.4e12*hdu.data.mean(axis=0), np.zeros_like(wavs),
             alpha=0.4, facecolor="k")
xlim(6500, 6700)
ylim(0.0, 0.3)
xlabel("Wavelength, Angstrom")
ylabel("Transmission")
legend()
gcf().set_size_inches(10, 8)

Tplot("FQ575N")
Tplot("F547M")
fill_between(wavs, 4e13*hdu.data.mean(axis=0), np.zeros_like(wavs),
             alpha=0.4, facecolor="k")
plot(*get_filter("F547M", 1, return_wavelength=True))
xlim(5000, 6000)
ylim(0.0, 0.3)
xlabel("Wavelength, Angstrom")
ylabel("Transmission")
legend()
gcf().set_size_inches(10, 8)

# F547M needs special treatment because it goes beyond the range of the red arm.   We simply extend Adal's spectrum down to 4900, using the mean value from the range 5450–5700

dwav, = np.diff(wavs[:2])

# Find indices corresponding to 5450–5700

i1 = np.argmin(np.abs(wavs - 5450.0))
i2 = np.argmin(np.abs(wavs - 5700.0))
i1, i2, wavs[102], wavs[404]

# Generate spatial profile of the mean continuum.

fill_values = np.mean(hdu.data[:, i1:i2], axis=1)

# Generate the extended wavelength range.

nx_extend = int((wavs[i1] - 4900)/dwav)
wavs_extend = wavs[i1] - (1. + np.arange(nx_extend)[::-1])*dwav
wavs_extend[0], wavs_extend[-1]

# Now stack the wavelength arrays together horizontally.

wavsx = np.hstack((wavs_extend, wavs[i1:]))

# The fill_values get stretched across the extended wavelength range. 

data_extend = np.ones(nx_extend)[None, :] * fill_values[:, None]
data_extend.shape

# And then they too are stacked horizontally with the existing data.

datax = np.hstack((data_extend, hdu.data[:, i1:]))
datax.shape, wavsx.shape

# Update the filters to use the extended wavelength range:

Tx = {f: get_interpolated_filter(f, wavsx) for f in wfc3_filters}
def Txplot(fname):
    plot(wavsx, Tx[fname], lw=3, alpha=0.6, label=fname)



# Now test it out:

Txplot("FQ575N")
Txplot("F547M")
mm = (wavsx > 6000.0) & (wavsx < 6200.0)
ph_datax = datax*wavsx[None, :]
meanspec = ph_datax.mean(axis=0)
for spec in ph_datax:
    plot(wavsx, 0.001*spec/spec[mm].mean(), 'k', lw=0.1, alpha=0.1)
floor = 5.e-4
plot(wavsx, 0.001*meanspec/meanspec[mm].mean(), 'y', alpha=0.7, lw=1.5)
xlim(5300, 6700)
ylim(floor, 1.0)
yscale("log")
xlabel("Wavelength, Angstrom")
ylabel("Transmission")
legend()
minorticks_on()
grid(ls='-', c='b', lw=0.3, alpha=0.3)
grid(ls='-', c='b', lw=0.3, alpha=0.1, which='minor')
gcf().set_size_inches(12, 8)
gcf().savefig("adal-red-slit6-lamFlam.png", dpi=600)

# Hurray – it worked!



# ### Now extract slit profiles for each of the filters.

filter_profiles = {f: np.sum(Tx[f][None, :]*datax, axis=1)*(wavsx[1]-wavsx[0])
                   for f in wfc3_filters}


xadal = 1.2*(xpos - 77.5)


# +
def Fplot(fname):
    plot(xadal, filter_profiles[fname], label=fname)
    
Fplot("F656N")
Fplot("F658N")
Fplot("F673N")
Fplot("F547M")
Fplot("FQ672N")
Fplot("FQ674N")
Fplot("FQ575N")
legend(ncol=2)
yscale('log')
xlabel('Slit position')
ylabel('Flux')
gcf().set_size_inches(8,6)

# -

# And plot some ratios too.

def Rplot(f1, f2, ax):
    ax.plot(xadal, filter_profiles[f1]/filter_profiles[f2])
    ax.set_ylabel("{} / {}".format(f1, f2))
    ax.grid()


# +
fig, axes = subplots(6, 1, sharex=True)
Rplot("FQ672N", "FQ674N", axes[0])
Rplot("F658N", "F656N", axes[1])
Rplot("FQ575N", "F658N", axes[2])
axes[2].set_ylim(0.0, 0.07)
Rplot("F673N", "F658N", axes[3])
axes[3].set_ylim(0.0, 0.6)
Rplot("F656N", "F547M", axes[4])
Rplot("FQ575N", "F547M", axes[5])

axes[-1].set_xlabel("Position, arcsec")
fig.set_size_inches(8,20)
fig.tight_layout()
# -

# ## Extract slit from WFC3 images

whdu = fits.open("full_F658N_north_pad.fits")[0]

wny, wnx = whdu.data.shape

import astropy.coordinates as coord
import astropy.units as u
from astropy import wcs
w = wcs.WCS(whdu.header)


# Define the position, orientation and width of Adal's slit

# +
def set_coord(coordstr):
    return coord.ICRSCoordinates(coordstr=coordstr, unit=(u.hour, u.deg))

slit6_center = set_coord("05:35:15.2 -05:23:53.1")
slit6_PA = 72.0
slit_width = 1.03
# -

# For each pixel in the WFC image, find the coordinate along the slit `xslit` and the coordinate across the slit `yslit`.

X, Y = np.meshgrid(np.arange(wnx), np.arange(wny))
RA, DEC = w.all_pix2world(X, Y, 0)
RA0, DEC0 = slit6_center.ra.deg, slit6_center.dec.deg
dRA = (RA - RA0)*3600*np.cos(np.radians(DEC0))
dDEC = (DEC - DEC0)*3600
COSPA, SINPA = np.cos(np.radians(slit6_PA)), np.sin(np.radians(slit6_PA))
xslit = dRA*SINPA + dDEC*COSPA
yslit = -dRA*COSPA + dDEC*SINPA

fits.PrimaryHDU(data=xslit, header=hdu.header).writeto("Adal_xslit6_north_pad.fits", clobber=True)
fits.PrimaryHDU(data=yslit, header=hdu.header).writeto("Adal_yslit6_north_pad.fits", clobber=True)


slitmask = np.abs(yslit) < 0.5*slit_width
print("Fraction of pixels in slit =", slitmask.sum()/(wny*wnx))

slitmask_p = np.abs(yslit+0.5) < 0.5*slit_width
slitmask_m = np.abs(yslit-0.5) < 0.5*slit_width

# +
plot(-xslit[slitmask_p]-1.5, whdu.data[slitmask_p], ",c", alpha=0.15)
plot(-xslit[slitmask_m]-1.5, whdu.data[slitmask_m], ",m", alpha=0.15)

plot(-xslit[slitmask]-1.5, whdu.data[slitmask], ",k", alpha=0.25)

plot(xadal, 1.6e13*filter_profiles["F658N"], "y", lw=3)
ylim(0.4, 30.0)
yscale("log")
xlim(-20.0, 90.0)
gcf().set_size_inches(12, 8)
# -

# Now smooth our image to match the ground-based resolution.

import scipy.ndimage as ni

pixscale = whdu.header["CDELT2"]*3600

clean_data = whdu.data[:, :]
mask = np.isfinite(clean_data)
clean_data[~mask] = np.median(clean_data[mask])

sdata = ni.gaussian_filter(clean_data, 1.0/pixscale)

# +
plot(-xslit[slitmask]-1.8, sdata[slitmask], ",k", alpha=0.25)

plot(xadal, 1.6e13*filter_profiles["F658N"], "or", lw=5)
ylim(0.4, 10.0)
#yscale("log")
xlim(-20.0, 90.0)
gcf().set_size_inches(12, 8)
# -

# Checking the slit against the images

sdata = Table.read("adal-slit6-EW.dat", format="ascii.tab")
s5data = Table.read("adal-slit5-EW.dat", format="ascii.tab")

sdata[np.isfinite(sdata["xim"])]

plot(sdata["Position"], sdata["xim"])

m5 = (s5data["xim"] > 0.0) & (s5data["xim"] < 60.0) 
m6 = (sdata["xim"] > -10.0) & (sdata["xim"] < 50.0) 

plot(sdata["R658"][m6], sdata["R575"][m6], 'o', alpha=0.3)
plot(s5data["R658"][m5], s5data["R575"][m5], 'or', alpha=0.3)
x = np.linspace(0.0, 8.0)
plot(x, 0.02 + 0.025*x)
xlim(0.0, None)
ylim(0.0, 0.3)

errorbar(s5data["xim"][m5], s5data["R575"][m5], yerr=s5data["R575e"][m5], fmt="ob")
errorbar(sdata["xim"][m6], sdata["R575"][m6], yerr=sdata["R575e"][m6], fmt="or")
ylim(0.0, 0.3)
gcf().set_size_inches(10, 6)
grid()

errorbar?

errorbar(s5data["Position"], s5data["R658"], yerr=s5data["R658e"], fmt='ob')
errorbar(sdata["Position"], sdata["R658"], yerr=sdata["R658e"], fmt='or')
ylim(0.0, None)
gcf().set_size_inches(10, 6)
grid()


# +
def Tm(T):
    "Maximum transmission of filter T"
    return T.max()

def Wj(wavs, T):
    "Rectangular width of filter T"
    return np.trapz(T, wavs)/Tm(T)

def Wtwid(wav0, wavs, T):
    """Find W-twiddle for a given line of wavelength wav0
    with respect to a filter transmission curve T(wavs)"""
    Ti = np.interp(wav0, wavs, T)
    # We are still missing the k_{j,i} term 
    return Tm(T)*Wj(wavs, T)/Ti


# +
import wfc3_utils

wav547, T547 = wfc3_utils.get_filter("F547M", return_wavelength=True)
wav575, T575 = wfc3_utils.get_filter("FQ575N", return_wavelength=True)
wav658, T658 = wfc3_utils.get_filter("F658N", return_wavelength=True)

def rqq(wav_line, wavs, TN, TW):
    r0 = Tm(TN)*Wj(wavs, TN)/(Tm(TW)*Wj(wavs, TW))
    q1 = (1./Wtwid(wav_line, wavs, TN) - 1./Wtwid(wav_line, wavs, TW))
    q2 = 1./(Wtwid(wav_line, wavs, TW) - Wtwid(wav_line, wavs, TN))
    return r0, q1, q2

r0, q1, q2 = rqq(5755, wav547, T575, T547)

def prelaunch_ratio(EW):
    return r0*(1.0 + q1*EW*(1.0 - q2*EW))
# -





r0, q1, q2 = rqq(6583, wav547, T658, T547)
plot(s5data["W6583"][m5], s5data["R658"][m5]/s5data["R547"][m5], 
     'o', alpha=0.6, label="Mesa-Delgado (2008) Slit 5")
plot(sdata["W6583"][m6], sdata["R658"][m6]/sdata["R547"][m6], 
     'or', alpha=0.6, label="Mesa-Delgado (2008) Slit 6")
x = np.linspace(0.0, 1200.0)
plot(x, prelaunch_ratio(x), label="Pre-launch calibration")
xlim(0.0, None)
ylim(0.0, None)
legend(loc="upper left")
gcf().set_size_inches(8, 6)
xlabel("EW(6583), Angstrom")
ylabel("R658 / (k R547)")
gcf().savefig("adal-w6583-calibration.pdf")



r0, q1, q2 = rqq(5755, wav547, T575, T547)
plot(s5data["W5755"][m5], s5data["R575"][m5]/s5data["R547"][m5]/0.93, 
     'o', alpha=0.6, label="Mesa-Delgado (2008) Slit 5")
plot(sdata["W5755"][m6], sdata["R575"][m6]/sdata["R547"][m6]/0.93, 
     'or', alpha=0.6, label="Mesa-Delgado (2008) Slit 6")
x = np.linspace(0.0, 25.0)
plot(x, prelaunch_ratio(x), label="Pre-launch calibration")
xlim(0.0, None)
ylim(0.0, 0.06)
legend(loc="lower right")
gcf().set_size_inches(8, 6)
xlabel("EW(5755), Angstrom")
ylabel("R575 / (k R547) ")
gcf().savefig("adal-w5755-calibration.pdf")

plot(s5data["xim"][m5], (1./600)*s5data["W6583"][m5], '-b')
plot(s5data["xim"][m5], s5data["R658"][m5]/s5data["R547"][m5], 'ob', alpha=0.6)
plot(sdata["xim"][m6], (1./600)*sdata["W6583"][m6], '-r')
plot(sdata["xim"][m6], sdata["R658"][m6]/sdata["R547"][m6], 'or', alpha=0.6)
xlim(None, None)
ylim(0.0, None)
gcf().set_size_inches(8, 6)

# ## Now try the ODH spectra

hdu = fits.open("ODell-Harris/S60.30.fits")[0]

nx, wav0, i0, dwav = [hdu.header[k] for k in ("NAXIS1", "CRVAL1", "CRPIX1", "CD1_1")]
wavs = wav0 + (np.arange(nx) - (i0 - 1))*dwav 

ny, dx = [hdu.header[k] for k in ("NAXIS2", "CD2_2")]
x0, j0 = 0.0, 0.0
xpos = x0 + (np.arange(ny) - (j0 - 1))*dx 

xpos

plot_profile(5754.6)
grid()

plot_profile(6547.9, trimred=5)
plot_profile(6583.4, width=18, trimblue=3)
plot_profile(6716.4, trimred=4)
plot_profile(6730.8, trimblue=4)
plot(wavs, hdu.data.mean(axis=0), lw=0.2)
yscale('log')
ylim(7.e-16, 8.e-13)
xlim(6450, 6900)
gcf().set_size_inches(12,8)

# +
wavs_547, f547m = get_filter("F547M", return_wavelength=True)
wavs_575, fq575n = get_filter("FQ575N", return_wavelength=True)
plot(wavs_547, 100*f547m, '--k', alpha=0.5, label='100 x T(F547M)')
fill_between(wavs_547, 100*f547m, np.zeros_like(wavs_547), alpha=0.05, zorder=-100, facecolor='k')
plot(wavs_575, 100*fq575n, ':k', alpha=0.8, label='100 x T(FQ575N)')
fill_between(wavs_575, 100*fq575n, np.zeros_like(wavs_575), alpha=0.1, zorder=-50, facecolor='k')

mm = (wavs > 6000.0) & (wavs < 6200.0)
ph_data = hdu.data*wavs[None, :]
meanspec = ph_data.mean(axis=0)
for spec in ph_data:
    plot(wavs, spec/spec[mm].mean(), 'k', lw=0.1, alpha=0.1)
floor = 2.e-1
plot(wavs, meanspec/meanspec[mm].mean(), 'y', alpha=0.7, lw=1.5, label="O'Dell & Harris (2010) mean spectrum")
xlim(3700, 7200)
ylim(floor, 100.0)
yscale("log")
xlabel("Wavelength, Angstrom")
ylabel("Transmission")
legend()
minorticks_on()
grid(ls='-', c='b', lw=0.3, alpha=0.3)
grid(ls='-', c='b', lw=0.3, alpha=0.1, which='minor')
gcf().set_size_inches(12, 8)
gcf().savefig("odh-S60-lamFlam.png", dpi=600)
# -

T["FQ575N"].shape, wavs.shape

T

annotate?


