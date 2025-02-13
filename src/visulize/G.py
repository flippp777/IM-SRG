import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import coo_array

data = np.transpose(np.loadtxt("output/G/Gamma_s0.0.txt"))
cols, rows, vals = data[0].astype(int), data[1].astype(int), np.sign(data[2])
data1 = np.transpose(np.loadtxt("output/G/Gamma_s0.3.txt"))
cols1, rows1, vals1 = data1[0].astype(int), data1[1].astype(int), np.sign(data1[2])
data2 = np.transpose(np.loadtxt("output/G/Gamma_s5.txt"))
cols2, rows2, vals2 = data2[0].astype(int), data2[1].astype(int), np.sign(data2[2])

N = int(140**2 *0.5)

def small_mat(n, cols, rows, vals):
    new_mat = np.zeros(shape=(n, n))
    i = 0
    len_data = len(cols)
    while (i < len_data):
        if cols[i] < n and rows[i] < n:
            new_mat[int(cols[i]), int(rows[i])] = vals[i]
        i = i+1
    return new_mat

# Gamma_small0 = small_mat(1000, cols, rows, vals)
# plt.figure()
# plt.imshow(Gamma_small0, cmap="seismic", vmax=5, vmin=-5)

# plt.figure()
# plt.spy(Gamma_small0, markersize=0.05)

# Gamma_small2 = small_mat(1000, cols2, rows2, vals2)
# plt.figure()
# plt.imshow(Gamma_small2, cmap="seismic", vmax=5, vmin=-5)

# plt.figure()
# plt.spy(Gamma_small2, markersize=0.05)


coo_Gamma = coo_array((vals, (rows, cols)))
coo_Gamma1 = coo_array((vals1, (rows1, cols1)))
coo_Gamma2 = coo_array((vals2, (rows2, cols2)))

# Gamma_mat2 = Gamma_mat[int(140*140/2):, int(140*140/2):]

plt.figure()
plt.spy(coo_Gamma, markersize=0.05)
plt.title("Gamma, s = 0")
plt.savefig("Gamma_s0.png", bbox_inches='tight')

plt.figure()
plt.spy(coo_Gamma1, markersize=0.05)
plt.title("Gamma, s = 0.3")
plt.savefig("Gamma_s0.3.png", bbox_inches='tight')

plt.figure()
plt.spy(coo_Gamma2, markersize=0.05)
plt.title("Gamma, s = 5")
plt.savefig("Gamma_s5.png", bbox_inches='tight')

# plt.figure()
# plt.imshow(coo_Gamma.toarray(), cmap="seismic", interpolation="none")
# plt.colorbar()
#
# plt.savefig("Gamma_mat.svg")

plt.show()
