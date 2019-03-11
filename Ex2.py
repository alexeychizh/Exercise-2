#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 20:36:46 2018

@author: alexchizh
"""
from numpy import *
import os, re
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

filepath = 'H2Ooutfiles' # Please feel free to change this to change between H2O and H2S
H2Xfiles = os.listdir(filepath)
XYZ = []
eq_geometry = zeros(3)
for file in H2Xfiles:
    rthetaE = re.findall(r'\d+\.\d+', file)
    f = list(open(os.path.join(filepath, file), 'r'))
    for line in f:
        if 'E(RHF)' in line:
            l = line.split()
            rthetaE.append(l[4])
    rthetaE_floats = []
    for i in rthetaE:
        rthetaE_floats.append(float(i))
    XYZ.append(rthetaE_floats)
    if rthetaE_floats[2] < eq_geometry[2]:
        for i in range(3):
            eq_geometry[i] = rthetaE_floats[i]
XYZ = asarray(XYZ)
X = XYZ[:, 0]
Y = XYZ[:, 1]
Z = XYZ[:, 2]
    
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_trisurf(X, Y, Z, cmap=cm.jet, linewidth=0)
fig.colorbar(surf)

ax.xaxis.set_major_locator(MaxNLocator(5))
ax.yaxis.set_major_locator(MaxNLocator(5))
ax.zaxis.set_major_locator(MaxNLocator(5))
ax.set_title('H2X energy landscape')
ax.set_xlabel('r / Å')
ax.set_ylabel('Theta / degrees')
ax.set_zlabel('Energy / Hartree')
fig.tight_layout()

plt.show()
print("Equilibrium Geometry:\n")
print("Lowest energy arrangement where r = %.2f Angstroms, theta = %.1f degrees\n" % (eq_geometry[0], eq_geometry[1]))
print("Energy of arrangement: %.3f Hartree" % (eq_geometry[2]))

""" vibrational frequencies """

x_eq = eq_geometry[0]
y_eq = eq_geometry[1]
E0 = eq_geometry[2]
mu = 1.660538921E-27

diff_theta = []
diff_r = []
diff_E1 = []
diff_E2 = []

for i in range(len(X)):
    if X[i] == x_eq:
        if abs((Y[i] - y_eq)) < 5:
            diff_theta.append((2*pi/360)*(Y[i] - y_eq))
            diff_E1.append(4.35974E-18*(Z[i] - E0))
    if Y[i] == y_eq:
        if abs((X[i] - x_eq)) < 0.2:
            diff_r.append(1E-10*(X[i] - x_eq))
            diff_E2.append(4.35974E-18*(Z[i] - E0))

p1 = polyfit(diff_theta, diff_E1, 2)
xfit1 = arange(-(2*pi/360)*5,(2*pi/360)*5, (2*pi/360)*0.1)
yfit1 = polyval(p1, xfit1)
plt.scatter(diff_theta, diff_E1, color='red', marker='+')
plt.plot(xfit1, yfit1)
plt.ylim(0, 4.35974E-21)
plt.title('Energy change on deviation from equilibrium Theta')
plt.xlabel('dTheta / degrees')
plt.ylabel('dEnergy/ Joules')
plt.show()

p2 = polyfit(diff_r, diff_E2, 2)
xfit2 = arange(-1E-10*0.2,1E-10*0.2, 1E-10*0.001)
yfit2 = polyval(p2, xfit2)
plt.scatter(diff_r, diff_E2, color='red', marker='+')
plt.plot(xfit2, yfit2)
plt.ylim(0, 1E-18)
plt.xlim(-1E-10*0.2,1E-10*0.2)
plt.title('Energy change on deviation from equilibrium r')
plt.xlabel('dr / m')
plt.ylabel('dEnergy/ Joules')
plt.show()

k_theta = 2*p1[0]
k_r = 2*p2[0]
vib_freq1 = (k_r/(2.0*mu))**0.5/(2*pi)
vib_freq2 = (k_theta/(0.5*mu*((1E-10*x_eq)**2)))**0.5/(2*pi)
vib_freq2 *= 1/2.9979E10
vib_freq1 *= 1/2.9979E10
print('Vibrational frequencies of normal modes:')
print('Bending mode: %.2f cm^-1' % (vib_freq2))
print('Symmetric stretching mode: %.2f cm^-1' % (vib_freq1))
