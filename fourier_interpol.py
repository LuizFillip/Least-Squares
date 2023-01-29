# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

x1 = np.array([100.,200.,230.,300.,400.,500.,600.,700.,800.,900.,1000, 1100, 1200])
y1 = np.array([21.87,38.34,695.47,762.60,695.47,351.41,180.49,102.88,64.57, 43.75,31.44, 23.66, 18.46])

def FourierEspelho(x1, y1):
    ndat = len(x1)
    n = ndat


    x = np.zeros(2 * ndat)
    y = np.zeros(2 * ndat)

    am = np.zeros(n)
    bm = np.zeros(n)
    rm = np.zeros(n)
    rm2 = np.zeros(n)
    fn = np.zeros(n)
    xmlat = np.zeros(n)
    sosf = np.zeros(n)

    for i in range(n):
        x[i] = x1[i]
        y[i] = y1[i]

    for i in range(ndat):
        j = ndat - i - 1
        y[ndat + i] = y[j]

    for m in range(n):
        nr = - n - 1
        sa = 0.0
        sb = 0.0
        for i in range(2*n):
            nr = nr + 1
            sa = sa + y[i] * np.cos(2 * np.pi * m * nr / (2 * n))
            sb = sb + y[i] * np.sin(2 * np.pi * m * nr / (2 * n))

        am[m] = (1.0 / (2*n)) * sa
        bm[m] = (1.0 / (2*n)) * sb
        rm[m] = np.sqrt(am[m]**2 + bm[m]**2) 
        rm2[m] = rm[m]**2

    for i in range(n):
        fn[i] = 1.0 * i / (2 * n * 0.5)

    f0 = 0
    m = 0

    resulty = []
    resultx = []

    step = 5
    for xpar in range(int(x1.min()), 
                      int(x1.max()) + step, 
                      step):

        if (xpar >= x1[-1]):
            xnr = - 1

        elif (xpar <= x1[0]):
            xnr = -ndat

        else:
            for i in range(ndat):
                if xpar == x[i]:
                    xnr = - ndat + i

                elif (xpar > x1[i]) and (xpar < x1[i + 1]):
                    xnr = (- ndat + i)*1

                    xnr = xnr + ((xpar - x1[i]) / (x1[i + 1] - x1[i]))


            st = 0
            for m in range(1, ndat):

                st += am[m]*np.cos(2*np.pi*xnr*0.5*m*fn[1]) + bm[m]*np.sin(2*np.pi*xnr*0.5*m*fn[1])

            resulty.append(am[0] + 2.0*st + am[-1]*np.cos(2*np.pi*xnr*0.5*n*fn[0]))
            resultx.append(xpar)
            
    return resultx, resulty 


fourier_x, fourier_y = FourierEspelho(x1, y1)

fig, ax = plt.subplots(figsize = (6, 8))
sns.set_style("whitegrid")

ax.plot(y1, x1, 'o', color = 'black', label = 'Dados')
ax.plot(fourier_y, fourier_x, label = 'Fourier')