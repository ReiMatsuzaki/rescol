import numpy as np
import matplotlib.pyplot as plt

for r0 in [50, 60, 70, 80, 90]:
    tmp = np.loadtxt("psi_{0}.dat".format(r0))
    plt.plot(tmp.T[0], tmp.T[1], label=str(r0))
plt.legend()
plt.savefig("psi.png")

plt.clf()

for r0 in [50, 60, 70, 80, 90]:
    tmp = np.loadtxt("psi_{0}.dat".format(r0))
    plt.plot(tmp.T[0], tmp.T[1], label=str(r0))
plt.legend()
plt.xlim(0, 5)
plt.savefig("psi_mini.png")

plt.clf()
for r0 in [50, 60, 70, 80, 90]:
    tmp = np.loadtxt("psi_{0}.dat".format(r0))
    plt.plot(tmp.T[0], tmp.T[1], label=str(r0))
plt.legend()
plt.xlim(40, 60)
plt.savefig("psi50.png")

