import heroespy
import matplotlib.pyplot as plt

files = heroespy.sas.get_all_files('ras')

ras = heroespy.sas.ras(files[30000])

fig, axes = plt.subplots(2, 1)
#fig = plt.subplots(nrows=2, sharex=True)
plt.subplot(211)
ras.plot()
plt.subplot(212)
plt.plot(ras.data[0,:])
for ax in axes:
    ax.set_anchor('W')

plt.show()