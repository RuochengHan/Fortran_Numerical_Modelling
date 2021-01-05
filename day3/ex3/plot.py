import numpy as np
import matplotlib.pyplot as plt

y = np.loadtxt(open('diff_result.dat', 'rb'))
x = np.array((list(range(1, 101))))
x = x / 100

axes = plt.gca()
axes.set_xlim([0,1])
axes.set_ylim([-2000,2000])

plt.plot(x, y)

plt.savefig("myplots.png")
