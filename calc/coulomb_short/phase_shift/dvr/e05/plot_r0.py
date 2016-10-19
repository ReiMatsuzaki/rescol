import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("r0_phase.dat")
plt.plot(data.T[0], data.T[1], "o")
print data.T[0]
print data.T[1]
plt.savefig("r0_phase.png")

