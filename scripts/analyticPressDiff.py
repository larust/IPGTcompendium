#!/usr/bin/python2
# -*- coding: utf-8 -*-
import CoolProp.CoolProp as cp
import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
#from matplotlib2tikz import save as tikz_save
# State variables
T=50+273.15
P0=1.0e5

# Propblem parameters
X=100.
b=10.
#k=9.869233e-13
k=1e-12
theta=0.5
P1=1.1e5
P2=1.0e5

mu=cp.PropsSI('V','T',T,'P',P0,'water')
rho=cp.PropsSI('D','T',T,'P',P0,'water')
beta=cp.DerivTermsU('IsothermalCompressibility',T,rho,'water',units='SI')

Deff=k/(theta*beta*mu)

N=100
n=np.array([[range(1,N+1)]])
m=np.array([[range(N)]])
x=np.transpose([[np.linspace(0,100,1000)]], (0,2,1))
t=np.transpose([[[2., 20., 200.]]],(2,0,1))

nMat=np.tile(n,(t.shape[0],x.shape[1],1))
mMat=np.tile(m,(t.shape[0],x.shape[1],1))
xMat=np.tile(x,(t.shape[0],1,n.shape[2]))
tMat=np.tile(t,(1,x.shape[1],n.shape[2]))

f=P1+(P2-P1)*xMat[:,:,0]/X+2/np.pi*np.sum((P2*np.cos(nMat*np.pi)-P1)/nMat*np.sin(nMat*np.pi*xMat/X)*np.exp(-Deff*np.power(nMat,2.)*np.pi**2*tMat/X**2),axis=2)+4*P0/np.pi*np.sum(1./(2*mMat+1)*np.sin((2*mMat+1)*np.pi*xMat/X)*np.exp(-Deff*np.power(2*mMat+1,2.)*np.pi**2*tMat/X**2),axis=2)

#ana=np.loadtxt('../Pressure diffusion BM/pressDiffAna.csv',delimiter=',',skiprows=1)
falc=np.loadtxt('../Pressure diffusion BM/pressDiffFALC.csv',delimiter=',',skiprows=1)
x=x.flatten()
p1 = plt.plot(x,f[0,:],'b-',x,f[1,:],'g-',x,f[2,:],'r-')
p2 = plt.plot(falc[:,0],falc[:,1],'bs',falc[:,0],falc[:,2],'gs',falc[:,0],falc[:,3],'rs')
plt.setp(p2, alpha=0.25)
plt.xlabel('Distance [m]')
plt.ylabel('Pressure [Pa]')
p1.extend(p2)
plt.legend(p1,['Analytic $t=2.00$ s','Analytic $t=20.0$ s','Analytic $t=200$ s','FALCON $t=2.00$ s','FALCON $t=20.0$ s','FALCON $t=200$ s'],loc='upper right')
#tikz_save('pressDiff.tikz',figureheight=r'\figureheight', figurewidth=r'\figurewidth',show_info=False)
plt.show()
