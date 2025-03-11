#!/usr/bin/env python
# -*- coding: utf-8 -*-
# gnana kotanja: metoda z vrtičnostjo v 2D
# TODO napaka: rešitev divergira

import numpy as np 
import matplotlib.pyplot as plt
from scipy import fft

plt.style.use('aljaz')

Re = 1 # Reynoldsovo število, brezdim

N = 50
L = 1 # dolžina stranice kvadrata
h = L/N # korak

cfl = h**2*Re/4 # CFL pogoj, max časovni korak
print(f"CFL = {cfl} s")
dt = 0.00005 # korak časa
Nmax = 10000 # max št korakov

tol = 1e-3 # toleranca konvergence

m=n=np.arange(1,1+N,1,dtype=int)
mm,nn=np.meshgrid(m,n) #indeksi

x=y=np.linspace(0,L,N)
X,Y=np.meshgrid(x,y) #mreža

# ZP
v = np.zeros((N,N)) # v_x
u = np.zeros((N,N)) # v_y
psi = np.zeros((N,N)) # tokovna funkcija, skalarna

v[0,:] = 1 # v_0 RP

zeta = np.zeros((N,N)) # vrtičnost, z-smer TODO boljši začetni približek!
for i in range(1,N-1):
    for j in range(1,N-1):
        zeta[i,j] = (u[j,i+1] - u[j,i-1] - (v[j+1,i] - v[j-1,i])) /2/h
u,v = u.T,v.T

# plt.imshow(zeta, cmap='viridis')
# plt.colorbar()
# plt.show()

fig,[ax1,ax2,ax,cax] = plt.subplots(1,4, figsize=(12,3), gridspec_kw={'width_ratios':[1,1,1,0.02]})
ax.set_aspect('equal')

# ITERACIJE
n=0
S_n=tol+1
conv=[]
sila=[]
while S_n > tol and n < Nmax:
    n+=1
    # 1. rešimo evolucijsko en za vrtičnost
    zeta_new = np.zeros((N,N))
    for i in range(1,N-1):
        for j in range(1,N-1):
            zeta_new[i,j] = zeta[i,j] - dt/2/h* (u[j,i+1]*zeta[j,i+1] - u[j,i-1]*zeta[j,i-1] + v[j+1,i]*zeta[j+1,i] - v[j-1,i]*zeta[j-1,i]) + \
                                        dt/h/h/Re* (psi[j,i+1] - 2*psi[j,i] + psi[j,i-1] + psi[j+1,i] - 2*psi[j,i] + psi[j-1,i])
    S_n = np.sum(np.abs(zeta-zeta_new)) #ustaljeno stanje, merilo
    conv.append(S_n)
    if S_n > 1e10: #divergira
        print('Divergenca: zmajšaj časovni korak!')
        break
    if n < 10: S_n = tol+1 #prvih nekaj korakov
    zeta = zeta_new

    # 2. rešimo Poissonovo enačbo za tokovno funkcijo - spektralno, alternative je MKR
    ZETA = fft.dstn(zeta, type=1) #transformiran prostor
    PSI = h*h/2 * ZETA / (np.cos(np.pi*mm/N)+np.cos(np.pi*nn/N)-2) #vektorsko
    psi = fft.idstn(PSI, type=1) #nazaj v fiz prostor

    # A=np.zeros((N*N,N*N)) #MKR - počasneje!
    # # DE
    # for j in range(1,N-1):
    #     for i in range(1,N-1):
    #         ij_index = j*(N) + i
    #         A[ij_index, ij_index] += 2*(1/h**2 + 1/h**2) #center
    #         # okolica
    #         A[ij_index, ij_index-1] += -1/h**2 #levo
    #         A[ij_index, ij_index+1] += -1/h**2 #desno
    #         A[ij_index, ij_index-N] += -1/h**2 #dol
    #         A[ij_index, ij_index+N] += -1/h**2 #gor      
    # # RP Dirichlet
    # for i in range(N):
    #     A[i, i] = 1 #spodaj
    #     A[-N+i, -N+i] = 1 #zgoraj
    # for j in range(N):
    #     A[j*N, j*N] = 1 #levo
    #     A[j*N+N-1, j*N+N-1] = 1 #desno
    # psi = np.linalg.solve(A, zeta.flatten()).reshape(N,N)
    
    # 3. hitrostno polje
    for i in range(1,N-1):
        for j in range(1,N-1):
            u[j,i] =  (psi[j+1,i] - psi[j-1,i])/2/h
            v[j,i] = -(psi[j,i+1] - psi[j,i-1])/2/h

    # divergenca hitrostnega polja
    # D = 0
    # for i in range(1,N-1):
    #     for j in range(1,N-1):
    #         D += abs(u[j,i+1] - u[j,i-1] + v[j+1,i] - v[j-1,i])/2/h /4

    # sila na pokrov
    f1 = h/Re * (zeta[0,0]/2 + np.sum(zeta[0,1:N-1]) + zeta[N-1,0]/2)
    f2 = h/Re * (zeta[N-1,0]/2 + np.sum(zeta[N-1,1:N-1]) + zeta[N-1,N-1]/2)
    f3 = h/Re * (zeta[0,0]/2 + np.sum(zeta[1:N-1,0]) + zeta[0,N-1]/2)
    f4 = h/Re * (zeta[N-1,0]/2 + np.sum(zeta[1:N-1,N-1]) + zeta[N-1,N-1]/2)    
    sila.append([f1,f2,f3,f4])

    if n % 100 == 0:
        V = np.sqrt(u**2 + v**2) #absolutna hitrost
        i=np.linspace(0,n,n)
        ax.clear()
        ax1.clear()
        ax2.clear()
        cax.clear()
        ax.grid(False)
        ax1.plot(i,sila,'-', label=['$f_1$','$f_2$','$f_3$','$f_4$'])
        ax1.legend(ncols=2, fontsize=10)
        ax2.plot(i,conv,'k-')
        # qv = ax.quiver(x,y, u,v, V, cmap='turbo')
        # fig.colorbar(qv, cax=cax, label='$v$ [m/s]')
        qv = ax.streamplot(x,y, u,v, color=V, cmap='turbo')
        fig.colorbar(qv.lines, cax=cax, label='$v$ [m/s]')
        fig.suptitle(f"$n = {n}$")
        fig.show()
        plt.pause(0.01)
        print(f"S_n = {S_n}")
print(f"Konvergenca po {n} korakih")
# print(f"Divergenca = {D}")
print(f"Ustaljeno stanje: {S_n}")
conv = np.asarray(conv)
sila = np.asarray(sila)
V = np.sqrt(u**2 + v**2) #absolutna hitrost
i=np.linspace(0,n,sila.shape[0])

# plt.contourf(X,Y,V, cmap='viridis')
ax.clear()
ax1.clear()
ax2.clear()
cax.clear()
ax.grid(False)
ax1.plot(i,sila,'-', label=['$f_1$','$f_2$','$f_3$','$f_4$'])
ax1.legend(ncols=2)
ax2.plot(i,conv,'k-')
# qv = ax.quiver(x,y, u,v, V, cmap='turbo')
# fig.colorbar(qv, cax=cax, label='$v$ [m/s]')
qv = ax.streamplot(x,y, u,v, color=V, cmap='turbo', broken_stremlines=False)
fig.colorbar(qv.lines, cax=cax, label='$v$ [m/s]')
fig.suptitle(f"$n = {n}$")
ax.set_xlabel('$x$ [m]')
ax.set_ylabel('$y$ [m]')
plt.tight_layout()
plt.show()

i=np.linspace(0,n,len(conv))
plt.plot(i,sila,'-', label=['$f_1$','$f_2$','$f_3$','$f_4$'])
# plt.xlabel('$t_n$ [s]')
plt.xlabel('$n$')
plt.ylabel('$F_n$')
plt.legend(ncols=2)
plt.show()

plt.plot(i,conv,'k-')
# plt.xlabel('$t_n$ [s]')
plt.xlabel('$n$')
plt.ylabel('$S_n$')
plt.show()