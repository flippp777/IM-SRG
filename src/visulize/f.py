import matplotlib.pyplot as plt
import numpy as np

fs001 = np.loadtxt("im_srg_euler_f_s0.txt")
fs018 = np.loadtxt("im_srg_euler_f_s0.0498.txt")
fs036 = np.loadtxt("im_srg_euler_f_s3.txt")

N = 140


def fill_mat(matrix_list):
    matrix = np.zeros(shape=(N, N))
    for element in matrix_list:
        matrix[int(element[0]), int(element[1])] = element[2]
    return matrix


fd001_mat = fill_mat(fs001)
fd018_mat = fill_mat(fs018)
fd036_mat = fill_mat(fs036)

plt.figure()
plt.title("f, s = 0")
plt.matshow(fd001_mat, fignum=0)
plt.colorbar()

plt.figure()
plt.title("f, s = 0.1")
plt.matshow(fd036_mat, fignum=0)
plt.colorbar()

plt.figure()
plt.title("s = 0.048")
plt.matshow(fd018_mat, fignum=0)

# plt.figure()
# plt.title("dfds, s = 0")
# plt.matshow(dfds_mat, fignum=0)
# plt.colorbar()

# plt.figure()
# plt.title("eta, s = 0")
# plt.matshow(eta_mat, fignum=0)
# plt.colorbar()

plt.show()
