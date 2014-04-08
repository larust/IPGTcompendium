#!/usr/bin/python2
# -*- coding: utf-8 -*-
import CoolProp.CoolProp as cp
import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
#from matplotlib2tikz import save as tikz_save
# State variables
T=50+273.15
P=1.0e5

# Propblem parameters
rout=100.
#k=9.869233e-13
k=1e-12
b=1.
beta=4.0e-10
theta=0.5
#g=9.80665
g=9.81
Qv=1e-3/60
rwell=1.
t=np.array([[2, 20, 200]]).T
# Physical properties determined from EoS
mu=cp.PropsSI('V','T',T,'P',P,'water')
rho=cp.PropsSI('D','T',T,'P',P,'water')
beta=cp.DerivTermsU('IsothermalCompressibility',T,rho,'water',units='SI')

# Other derived variables
Deff=k/(theta*beta*mu)
Tunit=k*b
Sunit=b*theta*beta
Kg=k*rho*g/mu
vinj=Qv/(2*np.pi*rwell*b)
Re=2*rwell*rho*vinj/mu

Trans = Tunit*rho*g/mu
Stor = Sunit*rho*g

r=np.linspace(rwell, rout, 1000)
rMat=np.tile(r,(len(t), 1))
tMat=np.tile(t,(1,len(r)))

u=np.power(rMat,2.)*Stor/(4*Trans*tMat)
W=sp.exp1(u)
s=Qv/(4*np.pi*Trans)*W
pres=P-s*rho*g

falc=np.loadtxt('../Theis problem BM/theis_falcon.csv',delimiter=',',skiprows=1)

p1 = plt.plot(r,pres[0,:],'b-',r,pres[1,:],'g-',r,pres[2,:],'r-')
p2 = plt.plot(falc[:,0],falc[:,1],'bs',falc[:,0],falc[:,2],'gs',falc[:,0],falc[:,3],'rs')
plt.setp(p2, alpha=0.25)
plt.xlabel('Radius [m]')
plt.ylabel('Pressure [Pa]')
p1.extend(p2)
plt.legend(p1,['Ana. $t=2.00$ s','Ana. $t=20.0$ s','Ana. $t=200$ s','FALCON $t=2.00$ s','FALCON $t=20.0$ s','FALCON $t=200$ s'],loc=4)
plt.ylim([9.2e4, 1e5])
#tikz_save('Theis.tikz',figureheight=r'\figureheight', figurewidth=r'\figurewidth',show_info=False)
plt.show()