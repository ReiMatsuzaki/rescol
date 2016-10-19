import numpy as np
import matplotlib.pyplot as plt

for num in [101, 201, 301]:
    tmp = np.loadtxt("psi_{0}.dat".format(num))
    plt.plot(tmp.T[0], tmp.T[1], label=str(num))
# plt.plot(tmp.T[0], tmp.T[2], "g")
# plt.xlim(44, 56)
# plt.ylim(-0.3, 0.7)
plt.legend()
plt.savefig("psi.png")

plt.clf()
for (num, c) in zip([101, 201, 301], ["r", "g", "b"]):
    tmp = np.loadtxt("psi_{0}.dat".format(num))
    plt.plot(tmp.T[0], tmp.T[1], c, label=str(num))
    plt.plot(tmp.T[0], tmp.T[2], c+"--", label=str(num))
plt.legend()
plt.xlim(70, 90)
plt.savefig("psi_mini.png")


