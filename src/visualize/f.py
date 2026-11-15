import matplotlib.pyplot as plt
import numpy as np

fs00 = np.loadtxt("output/f/im_srg_euler_f_s0.txt")
fs03 = np.loadtxt("output/f/f_s0.3.txt")
fs5 = np.loadtxt("output/f/f_s5.txt")

N = 140


def fill_mat(matrix_list):
    matrix = np.zeros(shape=(N, N))
    for element in matrix_list:
        matrix[int(element[0]), int(element[1])] = element[2]
    return matrix


f00_mat = fill_mat(fs00)
f03_mat = fill_mat(fs03)
f5_mat = fill_mat(fs5)

plt.figure()
plt.title("f, s = 0")
plt.matshow(f00_mat, fignum=0, cmap="seismic", vmin=-5, vmax=5)
plt.colorbar()
plt.savefig("f_s0.svg", bbox_inches='tight')

plt.figure()
plt.title("f, s = 0.3")
plt.matshow(f03_mat, fignum=0, cmap="seismic", vmin=-5, vmax=5)
plt.colorbar()
plt.savefig("f_s0.5.svg", bbox_inches='tight')

plt.figure()
plt.title("f, s = 5")
plt.matshow(f5_mat, fignum=0, cmap="seismic", vmin=-5, vmax=5)
plt.colorbar()
plt.savefig("f_s5.svg", bbox_inches='tight')

# plt.figure()
# plt.title("f, s = 0.3")
# plt.spy(f03_mat)

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

plt.show()
