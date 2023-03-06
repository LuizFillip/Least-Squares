# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 13:54:15 2021

@author: jussi
"""

import matplotlib.pyplot as plt 
import numpy as np


jmax, lmax = 20, 20

y1, y2 = 0, 1
x1, x2 = 0, 1

dx = (x2 - x1) / (lmax - 1)
dy = (y2 - y1) / (jmax - 1)

y = np.arange(y1, y2 + dy, dy)
x = np.arange(x1, x2 + dx, dx)

a = np.zeros((jmax,lmax))
b = np.zeros((jmax,lmax))
c = np.zeros((jmax,lmax))
d = np.zeros((jmax,lmax))
e = np.zeros((jmax,lmax))
f = np.zeros((jmax,lmax))
u_numerico = np.zeros((jmax,lmax)) 

for j in range(1, jmax):
    for l in range(1, lmax):
        

        a[j,l] = (1 / dy) * (1./dy) 
        b[j,l] = (1 / dy) * (1/dy) 
        c[j,l] = (1 / dx) * (1./dx) 
        d[j,l] = (1 / dx) * (1./dx) 
        e[j,l] = - 2*(1./dy)*(1./dy) - 2*(1./dx)*(1./dx) 
        f[j,l] = 0 
        u_numerico[j, l] = 0 

# COndições de contorno
u_numerico[:, 1] = np.sin(np.pi * x)
u_numerico[:, lmax - 1] = np.sin(np.pi * x)
u_numerico[1, :] = np.zeros(lmax)
u_numerico[jmax - 1, :] = np.zeros(lmax)

j1 = int(jmax / 2) + 1 
l1 = int(lmax / 2) + 1

f[j1, l1] = (1 / (jmax - 1)**2)  + (1  /  (lmax - 1)**2)


# Raio de Jacobi
jacobi_ratio = (np.cos(np.pi / jmax) + (dx /dy)**2 * np.cos(np.pi / lmax) ) / (1 + (dx / dy)**2)

maxits = 1000
eps = 0.00001
anormf= 0 
    
for j in range(2, jmax - 1):
    for l in range(2, lmax - 1):
        anormf = anormf + abs(f[j, l])
        
        
omega = 1

for n in range(1, maxits):
    anorm = 0
    jsw = 1
    for ipass in range(1, 2):
        lsw = jsw
        for j in range(2, jmax - 1):
            for l in np.arange(lsw + 1, lmax - 1, 2):
                resid = (a[j, l] * u_numerico[j + 1, 1] + b[j, l] * u_numerico[j - 1, 1] 
                         + c[j, l] * u_numerico[j, l + 1] + d[j, l] * u_numerico[j, l - 1] +
                         e[j, l]*u_numerico[j, l] - f[j, l])
                anorm = anorm + abs(resid)
                u_numerico[j, l] = u_numerico[j, l] - omega * resid / e[j, l]
            lsw - 3 - lsw
            
        jsw = 3 - jsw
        
        if (n == 1) and (ipass == 1):
            omega = 1 / (1 - 0.5 * jacobi_ratio**2)
        else:
            omega = 1 / (1 - 0.25 * jacobi_ratio**2 * omega)
 