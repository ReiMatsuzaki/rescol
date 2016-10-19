import numpy as np
import matplotlib.pyplot as plt

for num in range(40, 70):
   tmp = np.loadtxt("basis{0}.dat".format(num))
   plt.plot(tmp.T[0], tmp.T[1], "k")
   plt.plot(tmp.T[0], tmp.T[2], "g")
plt.xlim(44, 56)
plt.ylim(-0.3, 0.7)
plt.savefig("basis.png")

plt.clf()
tmp = np.loadtxt("basis60.dat".format(num))
plt.plot(tmp.T[0], tmp.T[1], "k")
plt.plot(tmp.T[0], tmp.T[2], "g")
plt.xlim(55, 65)
plt.ylim(-0.3, 0.7)
plt.savefig("basis60.png")


