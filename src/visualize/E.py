import matplotlib.pyplot as plt
import numpy as np

s1, E1, dEds1 = np.transpose(np.loadtxt("output/E/im_srg_euler_E_s00_s3_ds0.03.txt"))
srk, Erk = np.transpose(np.loadtxt("output/E/E_rk_h1e-06_s0.txt"))
s2, E2, dEds2 = np.transpose(np.loadtxt("output/E/im_srg_euler_E_s00_s3_ds0.001.txt"))


srk90, ErkN90 = np.transpose(np.loadtxt("output/E/E_rk_h1e-06_s0_N90.txt"))
srk40, ErkN40 = np.transpose(np.loadtxt("output/E/E_rk_h1e-06_s0_N40.txt"))

# plt.figure()
# plt.plot(s1, E1)
# plt.plot(s2, E2)
# plt.plot(srk, Erk)
# plt.title("Zero-Body Part of H for He4")
# plt.xlabel("s")
# plt.ylabel("E (MeV)")
# plt.grid(linestyle="--")
# plt.legend(["Euler: ds = 0.03", "Euler: ds = 0.001", "Runge Kutta: h=$1e-6$"])
# # plt.semilogx()
#
# plt.savefig("E_euler_rk.svg", bbox_inches='tight')

plt.figure()
plt.plot(srk, Erk)
plt.plot(srk90, ErkN90)
plt.plot(srk40, ErkN40)
plt.title("Zero-Body Part of H for He4, Runge-Kutta")
plt.xlabel("s")
plt.ylabel("E (MeV)")
plt.grid(linestyle="--")
plt.legend(["N=140", "N=90", "N=40"])
# plt.semilogx()

plt.savefig("E_rk_N.svg", bbox_inches='tight')

plt.show()
