from math import *
import numpy as np
from matplotlib import pyplot as plt

#Parameters

E = 1

Izz = 1

Iyy = 1

qx = 1

C1 = 1

C2 = 1

C3 = 1

C4 = 1

x1 =  0.153

x2 = 1.281

x3 = 2.681

xa = 0.28

P = 9

alpha = 0.5

F_y1 = 1

F_y2 = 1

F_y3 = 1

F_ya1 = 1

F_z1 = 1

F_z2 = 1

F_z3 = 1

F_za1 = 1

La = 2.771


locations = np.linspace(0,La)

#Deflection due to Mz in y direction

def Defl_Mz(x):
    v = -1/E/Izz*(-qx) + C1*x + C2

    if x > x1:
        v += -1/E/Izz*(F_y1/6*(x-x1)**(3))

    if x > x2-xa/2:
        v += -1/E/Izz*(F_ya1/6*(x-x2+xa/2)**(3))

    if x > x2:
        v += -1/E/Izz*(-F_y2/6*(x-x2)**(3))

    if x > x2+xa/2:
        v += -1/E/Izz*(-P*np.sin(alpha)/6*(x-x2-xa/2)**(3))

    if x > x3:
        v += -1/E/Izz*(F_y3/6*(x-x3)**(3))

    return v


valuesv = np.empty(50)
i = 0
for element in locations:
    valuesv[i] = Defl_Mz(element)
    i += 1

plt.plot(locations,valuesv)
plt.xlabel('Spanwise location on aileron [m]')
plt.ylabel('Deflection')
plt.title('Deflection in y direction due to bending moment Mz')
plt.show()

#Slope in y direction

def Slope_y(x):
    dvdx = -1/E/Izz*(-qx) + C1

    if x > x1:
        dvdx += -1/E/Izz*(F_y1/2*(x-x1)**(2))

    if x > x2-xa/2:
        dvdx += -1/E/Izz*(F_ya1/2*(x-x2+xa/2)**(2))

    if x > x2:
        dvdx += -1/E/Izz*(-F_y2/2*(x-x2)**(2))

    if x > x2+xa/2:
        dvdx += -1/E/Izz*(-P*np.sin(alpha)/2*(x-x2-xa/2)**(2))

    if x > x3:
        dvdx += -1/E/Izz*(F_y3/2*(x-x3)**(2))

    return dvdx

valuesdvdx = np.empty(50)
i = 0
for element in locations:
    valuesdvdx[i] = Slope_y(element)
    i += 1

plt.plot(locations,valuesdvdx)
plt.xlabel('Spanwise location on aileron [m]')
plt.ylabel('Slope')
plt.title('Slope in y direction due to bending moment Mz')
plt.show()


#Moment around z

def Moment_z(x):
    Mz = -qx

    if x > x1:
        Mz += F_y1*(x-x1)

    if x > x2-xa/2:
        Mz += F_ya1*(x-x2+xa/2)

    if x > x2:
        Mz += -F_y2*(x-x2)

    if x > x2+xa/2:
        Mz += -P*np.sin(alpha)*(x-x2-xa/2)

    if x > x3:
        Mz += F_y3*(x-x3)

    return Mz

valuesMz = np.empty(50)
i = 0
for element in locations:
    valuesMz[i] = Moment_z(element)
    i += 1

plt.plot(locations,valuesMz)
plt.xlabel('Spanwise location on aileron [m]')
plt.ylabel('Moment')
plt.title('Moment around z')
plt.show()

#Shear force in y

def Shear_y(x):
    Sy = -qx

    if x > x1:
        Sy += F_y1

    if x > x2-xa/2:
        Sy += F_ya1

    if x > x2:
        Sy += -F_y2

    if x > x2+xa/2:
        Sy += -P*np.sin(alpha)

    if x > x3:
        Sy += F_y3

    return Sy

valuesSy = np.empty(50)
i = 0
for element in locations:
    valuesSy[i] = Shear_y(element)
    i += 1

plt.plot(locations,valuesSy)
plt.xlabel('Spanwise location on aileron [m]')
plt.ylabel('Shear Force')
plt.title('Shear force in y direction')
plt.show()


#Deflection due to My in z direction

def Defl_My(x):
    w = C3*x + C4

    if x > x1:
        w += -1/E/Izz*(F_z1/6*(x-x1)**(3))

    if x > x2-xa/2:
        w += -1/E/Izz*(F_za1/6*(x-x2+xa/2)**(3))

    if x > x2:
        w += -1/E/Izz*(F_z2/6*(x-x2)**(3))

    if x > x2+xa/2:
        w += -1/E/Izz*(-P*np.cos(alpha)/6*(x-x2-xa/2)**(3))

    if x > x3:
        w += -1/E/Izz*(F_z3/6*(x-x3)**(3))

    return w

valuesDefl_My = np.empty(50)
i = 0
for element in locations:
    valuesDefl_My[i] = Defl_My(element)
    i += 1
    
plt.plot(locations,valuesDefl_My)
plt.xlabel('Spanwise location on aileron [m]')
plt.ylabel('Deflection')
plt.title('Deflection in z direction due to bending moment My')
plt.show()

#slope in z direction

def Slope_z(x):
    dwdx = C3

    if x > x1:
        dwdx += -1/E/Iyy*(F_z1/2*(x-x1)**(2))

    if x > x2-xa/2:
        dwdx = -1/E/Iyy*(F_za1/2*(x-x2+xa/2)**(2))

    if x > x2:
        dwdx += -1/E/Iyy*(F_z2/2*(x-x2)**(2))

    if x > x2+xa/2:
        dwdx += -1/E/Iyy*(-P*np.cos(alpha)/2*(x-x2-xa/2)**(2))

    if x > x3:
        dwdx += -1/E/Iyy*(F_z3/2*(x-x3)**(2))

    return dwdx

valuesdwdx = np.empty(50)
i = 0
for element in locations:
    valuesdwdx[i] = Slope_z(element)
    i += 1

plt.plot(locations,valuesdwdx)
plt.xlabel('Spanwise location on aileron [m]')
plt.ylabel('Slope')
plt.title('Slope in z direction due to bending moment My')
plt.show()

#Moment around y

def Momenty(x):
    My = 0

    if x > x1:
        My += F_z1*(x-x1)

    if x > x2-xa/2:
        My += F_za1*(x-x2+xa/2)

    if x > x2:
        My += F_z2*(x-x2)

    if x > x2+xa/2:
        My += -P*np.cos(alpha)*(x-x2-xa/2)

    if x > x3:
        My += F_z3*(x-x3)

    return My

valuesMy = np.empty(50)
i = 0
for element in locations:
    valuesMy[i] = Momenty(element)
    i += 1

plt.plot(locations,valuesMy)
plt.xlabel('Spanwise location on aileron [m]')
plt.ylabel('Moment')
plt.title('Moment around y')
plt.show()

#shear force in z direction

def Shear_z(x):
    Sz = 0

    if x > x1:
        Sz += F_z1

    if x > x2-xa/2:
        Sz += F_za1

    if x > x2:
        Sz += F_z2

    if x > x2+xa/2:
        Sz += -P*np.cos(alpha)

    if x > x3:
        Sz += F_z3

    return Sz

valuesSz = np.empty(50)
i = 0
for element in locations:
    valuesSz[i] = Shear_z(element)
    i += 1

plt.plot(locations,valuesSz)
plt.xlabel('Spanwise location on aileron [m]')
plt.ylabel('Shear Force')
plt.title('Shear force in y direction')
plt.show()
