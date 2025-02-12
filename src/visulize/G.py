import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import coo_array

data = np.transpose(np.loadtxt("output/G/Gamma_s0.0.txt"))
cols, rows, vals = data[0].astype(int), data[1].astype(int), np.sign(data[2])
N = int(140**2 *0.5)

coo_Gamma = coo_array((vals, (rows, cols)))


def fill_mat(matrix_list):
    matrix = np.empty(shape=(N, N))
    matrix[:] = np.nan
    for element in matrix_list:
        if int(element[0]) < N and int(element[1]) < N:
            matrix[int(element[0]), int(element[1])] = element[2]
    return matrix


# Gamma_mat2 = Gamma_mat[int(140*140/2):, int(140*140/2):]

plt.figure()
plt.spy(coo_Gamma, markersize=0.05)
plt.title("Gamma, s = 0")

# plt.figure()
# plt.imshow(coo_Gamma.toarray(), cmap="seismic", interpolation="none")
# plt.colorbar()
#
# plt.savefig("Gamma_mat.svg")

plt.show()
