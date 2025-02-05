import numpy as np
import matplotlib.pyplot as plt

fs001 = np.loadtxt("output/im_srg_euler_f_s0.txt")
fs018 = np.loadtxt("output/im_srg_euler_f_s0.018.txt")
fs036 = np.loadtxt("output/im_srg_euler_f_s0.036.txt")
fs063 = np.loadtxt("output/im_srg_euler_f_s0.063.txt")

N = 140
def fill_mat(matrix_list):
    matrix = np.zeros(shape=(N, N))
    for element in matrix_list:
        matrix[int(element[0]), int(element[1])] = element[2]
    return matrix

fd001_mat = fill_mat(fs001)
fd018_mat = fill_mat(fs018)
fd036_mat = fill_mat(fs036)
fd063_mat = fill_mat(fs063)

plt.figure()
plt.title("s = 0")
plt.matshow(fd001_mat, fignum = 0)

plt.figure()
plt.title("s = 0.018")
plt.matshow(fd018_mat, fignum = 0)

plt.figure()
plt.title("s = 0.036")
plt.matshow(fd036_mat, fignum = 0)

plt.show()
