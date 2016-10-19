import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("ene_phase.dat")
plt.plot(data.T[0], data.T[1], "o")
plt.savefig("phase.png")

