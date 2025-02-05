import matplotlib as plt
import numpy as np

s, E, dEds = np.loadtxt("output/im_srt_euler_E.txt")

plt.figure()
plt.plot(s, E, "o-")
plt.xlable("s")
plt.ylabel("E (MeV)")
plt.gird(linestyle="--")

plt.figure()
plt.plot(s, dEds, "o-")
plt.xlable("s")
plt.ylabel("dEds (MeV)")
plt.gird(linestyle="--")
