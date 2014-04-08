#!/usr/bin/python2
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib2tikz import save as tikz_save
rho=7850.
E=2e11
G=.3e11
nu=0.225
L=1.
h=.01
q=1e3
g=9.80665

F=q*L
D=E*np.power(h,3)/(12*(1-np.power(nu,2)))
I=D/E

bf=np.loadtxt('../Beam deflection BM/M1_beam_deflection.csv',delimiter=',',skiprows=1)

x=np.linspace(0,1.,1000)
M=.5*q*(L-x)**2
w=-q*x**2/D*(x**2/24. - L*x/6. +L**2/4)
sigmaMax=-6*M/h**2

fig1, ax1 = plt.subplots()
ax1.plot(x,w*1e3,'r-')
#ax2.plot(x,sigmaMax/q/1e3,'g-')
ax1.set_xlabel('Fractional position along plate of length L')
ax1.set_ylabel('Displacement [mm]')
#ax2.set_ylabel('Fiber stress, relative to loading [$\cdot10^3$]')
tikz_save('beamDisp.tikz',figureheight=r'\figureheight', figurewidth=r'\figurewidth',show_info=False)

fig2, ax2 = plt.subplots()
ax2.plot(x,sigmaMax/q/1e3)
ax2.set_ylabel('Fiber stress, relative to loading [$\cdot10^3$]')
ax2.set_xlabel('Fractional position along plate of length L')
#tikz_save('beamStress.tikz',figureheight=r'\figureheight', figurewidth=r'\figurewidth',show_info=False)
plt.show()
