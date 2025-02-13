import matplotlib.pyplot as plt
import numpy as np

etas00 = np.loadtxt("output/eta/eta1_s0.txt")
etas03 = np.loadtxt("output/eta/eta1_s0.3.txt")
etas5 = np.loadtxt("output/eta/eta1_s5.txt")

N = 140


def fill_mat(matrix_list):
    matrix = np.zeros(shape=(N, N))
    for element in matrix_list:
        matrix[int(element[0]), int(element[1])] = element[2]
    return matrix


feta00_mat = fill_mat(etas00)
eta03_mat = fill_mat(etas03)
eta5_mat = fill_mat(etas5)

plt.figure()
plt.title("eta, s = 0")
plt.matshow(feta00_mat, fignum=0, cmap="seismic", vmin=-5, vmax=5)
plt.colorbar()
plt.savefig("eta_s0.svg", bbox_inches='tight')

plt.figure()
plt.title("eta, s = 0.3")
plt.matshow(eta03_mat, fignum=0, cmap="seismic", vmin=-5, vmax=5)
plt.colorbar()
plt.savefig("eta_s0.5.svg", bbox_inches='tight')

plt.figure()
plt.title("eta, s = 5")
plt.matshow(eta5_mat, fignum=0, cmap="seismic", vmin=-5, vmax=5)
plt.colorbar()
plt.savefig("eta_s5.svg", bbox_inches='tight')


plt.show()
