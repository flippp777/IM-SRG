import matplotlib.pyplot as plt
import numpy as np

s, E, dEds = np.transpose(np.loadtxt("output/im_srg_euler_E.txt"))

plt.figure()
plt.plot(s, E, "o-")
plt.xlim(right = 0.04)
plt.xlabel("s")
plt.ylim(bottom = -30)
plt.ylabel("E (MeV)")
plt.grid(linestyle="--")

plt.figure()
plt.plot(s, dEds, "o-")
plt.xlim(right = 0.04)
plt.ylim(bottom = -200)
plt.xlabel("s")
plt.ylabel("dEds (MeV)")
plt.grid(linestyle="--")

plt.show()
