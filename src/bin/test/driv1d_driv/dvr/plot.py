import numpy as np
import matplotlib.pyplot as plt

for (num, c) in zip([101, 201, 301], ["r", "g", "b"]):
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
plt.xlim(60, 80)
plt.savefig("psi_mini.png")

plt.clf()
for (num, c) in zip([101, 201, 301], ["r", "g", "b"]):
    tmp = np.loadtxt("psi_{0}.dat".format(num))
    plt.plot(tmp.T[0], tmp.T[1], c, label=str(num))
    plt.plot(tmp.T[0], tmp.T[2], c+"--", label=str(num))
plt.legend()
plt.xlim(60, 62)
plt.savefig("psi_super_mini.png")

plt.clf()
for (num, c) in zip([101, 201, 301], ["r", "g", "b"]):
    tmp = np.loadtxt("psi_{0}.dat".format(num))
    plt.plot(tmp.T[0], tmp.T[1], c, label=str(num))
    plt.plot(tmp.T[0], tmp.T[2], c+"--", label=str(num))
plt.legend()
#tmp = np.loadtxt("psi_101.dat")
#plt.plot(tmp.T[0], tmp.T[1])
#plt.plot(tmp.T[0], tmp.T[2])
plt.xlim(0, 5)
plt.savefig("psi_near0.png")


