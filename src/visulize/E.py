import matplotlib.pyplot as plt
import numpy as np

# s, E, dEds = np.transpose(np.loadtxt("output/im_srg_euler_E_s00_s3_ds0.03.txt"))
# s, E = np.transpose(np.loadtxt("E_runge-kutta.txt"))
s, E, dEds = np.transpose(np.loadtxt("im_srg_euler_E_s00_s3_ds0.001.txt"))


# i = 0
# new_s = []
# while (i<len(s)):
#     new_s.append(i*0.02)
#     i += 1
#
# s = new_s

plt.figure()
plt.plot(s, E, "o-")
# plt.xlim(right = 0.04)
plt.xlabel("s")
# plt.ylim(bottom = -30)
plt.ylabel("E (MeV)")
plt.grid(linestyle="--")
# plt.semilogx()

# plt.figure()
# plt.plot(s, dEds, "o-")
# # plt.ylim(bottom = -200)
# plt.xlabel("s")
# plt.ylabel("dEds (MeV)")
# plt.grid(linestyle="--")

plt.show()
