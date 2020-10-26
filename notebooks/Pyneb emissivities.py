# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %matplotlib inline
import numpy as np
import pyneb as pn
from matplotlib import pyplot as plt
import seaborn as sns


oiir = pn.RecAtom("O", 2)

oiiic = pn.Atom("O", 3)

oiiic.getTransition(5007)

pn.atomicData.getDataFile('O2', 'rec')

pn.atomicData.getAllAtoms(coll=False, rec=True)

[_ for _ in oiir.labels if _.startswith("43")]

oiir.getTransition('4649.13')

oiir.getEmissivity(tem=np.linspace(5000, 10000, 6), den=1e4, label='4650.84')

temperatures = np.linspace(1000.0, 15000.0, 200)

# +
V1_mult = ('4638.86', '4641.81', '4649.13', '4650.84', '4661.63', '4673.73', '4676.23', '4696.35')
def oii_multiplet_sum_emissivity(tem, den, multiplet=V1_mult):
    result = None
    for label in multiplet:
        em = oiir.getEmissivity(tem=tem, den=den, label=label)
        if result is None:
            result = em
        else:
            result += em
    return result


def normalize(tem, em, T0=8000.0):
    # Normalize to T = T0
    em0 = np.interp(T0, tem, em)
    return em / em0
    
 
# -

V15_mult = "4590.97", #"4596.18", "4595.96"
V3d_4f = "4303.82", "4609.44", "4087.15", "4089.29"

oii_multiplet_sum_emissivity([1000, 3000, 5000, 8000, 1e4], 1e4, multiplet=V15_mult)

sns.set_context("poster")
sns.set_color_codes("bright")

# +
fig, ax = plt.subplots(figsize=(8, 8))
density = 1e4

em = oii_multiplet_sum_emissivity(temperatures, density)
em = normalize(temperatures, em)
ax.plot(temperatures, em, label="O II V1", ls="-", lw=6)

em = oii_multiplet_sum_emissivity(temperatures, density, multiplet=V15_mult)
em = normalize(temperatures, em)
ax.plot(temperatures, em, label="O II V15", ls="--", lw=6)

em = oii_multiplet_sum_emissivity(temperatures, density, multiplet=V3d_4f)
em = normalize(temperatures, em)
ax.plot(temperatures, em, label="O II 3d-4f", ls=":", lw=6)

em = oiiic.getEmissivity(tem=temperatures, den=density, wave=5007)
em = normalize(temperatures, em)
ax.plot(temperatures, em, label="[O III] 5007", ls="-", lw=2.5)

em = oiiic.getEmissivity(tem=temperatures, den=density, wave=4363)
em = normalize(temperatures, em)
ax.plot(temperatures, em, label="[O III] 4363", ls="-.", lw=2.5)

ax.legend(fontsize="x-small")
ax.set(
    xscale="linear",
    yscale="linear",
    ylim=[0, 3],
    xlim=[0, None],
    xlabel="Temperature, K",
    ylabel="Relative emissivity",
)
sns.despine()
fig.tight_layout()
fig.savefig("oii-emissivity-vs-t.pdf")

# +
# oiiic.getEmissivity?
