import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

data = np.loadtxt(open('result.dat', 'rb'))

fig = plt.figure(1)
cmap = mpl.colors.LinearSegmentedColormap.from_list('my_colormap',
                                           ['blue', 'yellow', 'red'],
                                           256)
img = plt.imshow(data,interpolation='nearest',
                 cmap = cmap,
                 origin='lower')
plt.colorbar(img, cmap=cmap)
fig.savefig("image.png")
