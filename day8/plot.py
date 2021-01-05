import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import f90nml
import os

parser = argparse.ArgumentParser()
parser.add_argument('-d', type=str, default=None, help='-d [.dat file name]')
args = parser.parse_args()

data = np.loadtxt(open(args.d, 'rb'))

fig = plt.figure(1)
cmap = mpl.colors.LinearSegmentedColormap.from_list('my_colormap',
                                           ['blue', 'red'],
                                           256)
img = plt.imshow(data, interpolation='nearest',
                 cmap = cmap,
                 origin='lower')
plt.colorbar(img, cmap=cmap)

nml = f90nml.read('params.dat').get("inputs")
Pr = str(nml.get("pr"))
t = str(nml.get("total_time"))
Ra = str(nml.get("ra"))

path = "./Pr-" + Pr + "_time-" + t + "_Ra-" + Ra
try:
    os.stat(path)
except:
    os.mkdir(path)

filename = args.d.replace(".dat", ".png")
filename = os.path.join(path, filename)

fig.savefig(filename)
