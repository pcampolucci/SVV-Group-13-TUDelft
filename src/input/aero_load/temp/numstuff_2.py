from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

loads = np.genfromtxt("load_A380.dat", delimiter=",")
testloads = np.ones([81, 41])

zset = []
xset = []
ca = 0.547
la = 2.771
for i in range(1, 82):
    thetaz = (i - 1) * np.pi / 81
    thetaz1 = i * np.pi / 81
    z = -0.5 * (0.5 * ca * (1 - np.cos(thetaz)) + 0.5 * ca * (1 - np.cos(thetaz1)))
    zset.append(z)

for i in range(1, 42):
    thetaz = (i - 1) * np.pi / 41
    thetaz1 = i * np.pi / 41
    x = 0.5 * (0.5 * la * (1 - np.cos(thetaz)) + 0.5 * la * (1 - np.cos(thetaz1)))
    xset.append(x)

zstuff = np.array(zset)
xstuff = np.array(xset)


"---------- Linear spline interpolation -----------"

h_z = np.diff(zstuff)
h_aload = np.diff(loads.T[0])

coeffs = np.zeros((80, 2))

coeffs.T[0] = loads.T[0][:-1]
coeffs.T[1] = h_aload / h_z

zvalue = zstuff[0]
for i in range(len(zstuff) - 1):
    loadsety = []
    realz = []
    while zvalue > zstuff[i + 1]:
        loadiny = coeffs.T[0][i] + coeffs.T[1][i] * (zvalue - zstuff[i])
        loadsety.append(loadiny)
        realz.append(zvalue)
        zvalue -= 0.00005

    plt.plot(realz, loadsety)

plt.scatter(zstuff, loads.T[0])
plt.show()
