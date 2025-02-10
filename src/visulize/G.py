import matplotlib.pyplot as plt
import numpy as np

Gamma = np.loadtxt("output/G/Gamma_s0.0.txt")
N = int(140**2*0.75)


def fill_mat(matrix_list):
    matrix = np.empty(shape=(N, N))
    matrix[:] = np.nan
    for element in matrix_list:
        if int(element[0]) < N and int(element[1]) < N:
            matrix[int(element[0]), int(element[1])] = element[2]
    return matrix

Gamma_mat = fill_mat(Gamma)
# Gamma_mat2 = Gamma_mat[int(140*140/2):, int(140*140/2):]

plt.figure()
plt.title("Gamma, s = 0")
plt.matshow(Gamma_mat, fignum=0)
plt.colorbar()


# plt.figure()
# plt.title("s = 0.048")
# plt.matshow(fd018_mat, fignum=0)

# plt.figure()
# plt.title("dfds, s = 0")
# plt.matshow(dfds_mat, fignum=0)
# plt.colorbar()

# plt.figure()
# plt.title("eta, s = 0")
# plt.matshow(eta_mat, fignum=0)
# plt.colorbar()
plt.savefig("Gamma_mat.pdf")

plt.show()
