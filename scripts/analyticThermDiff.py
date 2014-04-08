#!/usr/bin/python2
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
#from matplotlib2tikz import save as tikz_save
lf=0.68
lr=2.955
cpf=4186.
cpr=920.
rhof=1000.
rhor=2500.
T0=0.
Tinj=1.
X=10.
t1=1e6
t2=1e7
theta=0.2

leff=lf*theta+lr*(1-theta)
ceff=rhof*cpf*theta+rhor*cpr*(1-theta)
Deff=leff/(rhof*cpf)
R=ceff/(theta*rhof*cpf)

Deff=Deff/(R*theta)
x=np.linspace(0,X,1000)
T=erf(x/2./np.sqrt(Deff*t1))
T=np.vstack((T,erf(x/2./np.sqrt(Deff*t2))))

bf=np.loadtxt('../Heat CDE BM/HeatDiffusionBM-FALCON.csv',delimiter=',',skiprows=1)
xF=bf[:,0]
TF1=bf[:,1]
TF2=bf[:,2]

p1=plt.plot(x,T[0,:],'r-',x,T[1,:],'g-')
p2=plt.plot(xF,TF1,'rs',xF,TF2,'gs')
plt.setp(p2,alpha=0.25)
plt.ylim([0., 1.])
plt.xlabel('Distance [m]')
plt.ylabel('Temperature [$^\circ$C]')
p1.extend(p2)
plt.legend(p1,['Analytic, $t=10^6$ s', 'Analytic, $t=10^7$ s','FALCON, $t=10^6$ s', 'FALCON, $t=10^7$ s'],loc=4)
#tikz_save('heatDiff.tikz',figureheight=r'\figureheight', figurewidth=r'\figurewidth',show_info=False)
plt.show()
