import matplotlib.pyplot as plt
import numpy as np

s1, E1, dEds1 = np.transpose(np.loadtxt("output/E/im_srg_euler_E_s00_s3_ds0.03.txt"))
srk, Erk = np.transpose(np.loadtxt("output/E/E_runge-kutta.txt"))
s2, E2, dEds2 = np.transpose(np.loadtxt("output/E/im_srg_euler_E_s00_s3_ds0.001.txt"))


# i = 0
# new_s = []
# while (i<len(s)):
#     new_s.append(i*0.02)
#     i += 1
#
# s = new_s

plt.figure()
plt.plot(s1, E1)
plt.plot(s2, E2)
plt.plot(srk, Erk)
plt.xlabel("s")
plt.ylabel("E (MeV)")
plt.grid(linestyle="--")
plt.legend(["ds = 0.03", "ds = 0.001", "rk"])
plt.semilogx()

# plt.figure()
# plt.plot(s, dEds, "o-")
# # plt.ylim(bottom = -200)
# plt.xlabel("s")
# plt.ylabel("dEds (MeV)")
# plt.grid(linestyle="--")

plt.show()
