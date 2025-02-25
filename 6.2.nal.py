#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 6. naloga: PDE lastne rešitve
# homogena polkrožna opna

import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer
from numba import njit
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh
import matplotlib.cm as cm
import matplotlib.colors as mcolors

plt.style.use('aljaz')
np.set_printoptions(linewidth=200, precision=3, suppress=True)

@njit(cache=True)
def inv_potenčna(A, sigma, tol=1e-9,Nmax=10000):
    """inverzna potenčna metoda za iskanje lastne vrednosti matrike A blizu sigma"""
    A_inv=np.linalg.inv(A-sigma*np.empty_like(A))
    x=np.random.rand(A.shape[0])
    l=0
    l_old=1
    i=0
    while i<Nmax and abs(l-l_old)>tol:
        l_old=l
        y=A_inv@x
        x=y/np.linalg.norm(y,2)
        l=x.T@A@x
        i+=1
    print('lambda =', round(l,3), ', iteracij = ', i, ', tol = ', abs(l-l_old))
    return l, x
def homog(N):
    N_r=N #št. računskih točk v radialni smeri
    N_phi=N #št. računskih točk v obodni smeri
    R=1 #polmer
    PHI = np.pi
    h_r=R/(N_r-1)
    h_phi=PHI/(N_phi-1)

    phi=np.linspace(0,PHI,N_phi)
    r=np.linspace(0,R,N_r)
    r,phi=np.meshgrid(r,phi)
    x=r*np.cos(phi)
    y=r*np.sin(phi)

    # matrika A
    A = np.zeros((N_r*N_phi,N_r*N_phi))
    for i in range(N_r):
        for j in range(N_phi):
            A[i+j*N_phi, i+j*N_phi] = 2*(1/h_r**2 + 1/(h_r**2*h_phi**2)*1/(i+1)**2)  # Main diagonal i,j
            if i > 0:
                A[i+j*N_phi, i-1+j*N_phi] = 1/h_r**2 + 1/(2*(i+1)*h_r**2)  # i-1,j
            if j > 0:
                A[i+j*N_phi, i+(j-1)*N_phi] = 1/(h_r**2*h_phi**2*(i+1)**2)  # i,j-1
            if i < N_phi-1:
                A[i+j*N_phi, i+1+j*N_phi] = 1/h_r**2 + 1/(2*(i+1)*h_r**2)  # i+1,j
            if j < N_phi-1:
                A[i+j*N_phi, i+(j+1)*N_phi] = 1/(h_r**2*h_phi**2*(i+1)**2)  # i,j+1
    # plt.imshow(A)
    # plt.show()

    # problem last vrednosti
    # w1, x1 = np.linalg.eig(A)
    t1=timer()
    # w, X = np.linalg.eigh(A) #sym -> hitreje!
    w, X = eigh(A) #scipy
    t_eigh=timer()-t1
    print('t eigh = ', round(t_eigh,3), ' s')
    t1=timer()
    ws, Xs = eigsh(A, k=24, which='LM', sigma=0) #sparse shift invert mode!
    t_eigsh=timer()-t1
    print('t eigsh = ', round(t_eigsh,3), ' s')
    # print(np.sum(np.abs(w1-w)), np.sum(np.abs((x1-x)))) #ok za nizke lastne vrednosti
    t1=timer()
    wp,xp=inv_potenčna(A, 0)
    tp=timer()-t1
    return (w,ws,wp), (X,Xs,xp), (t_eigh,t_eigsh,tp), (x,y)
def nehomog(rho):
    N=rho.shape[0]
    N_r=N #št. računskih točk v radialni smeri
    N_phi=N #št. računskih točk v obodni smeri
    R=1 #polmer
    PHI = np.pi
    h_r=R/(N_r-1)
    h_phi=PHI/(N_phi-1)

    phi=np.linspace(0,PHI,N_phi)
    r=np.linspace(0,R,N_r)
    r,phi=np.meshgrid(r,phi)
    x=r*np.cos(phi)
    y=r*np.sin(phi)

    # matrika A
    A = np.zeros((N_r*N_phi,N_r*N_phi))
    B=np.zeros_like(A)
    for i in range(N_r):
        for j in range(N_phi):
            B[i+j*N, i+j*N] = rho[i,j] #desna stran
            A[i+j*N_phi, i+j*N_phi] = 2*(1/h_r**2 + 1/(h_r**2*h_phi**2)*1/(i+1)**2)  # Main diagonal i,j
            if i > 0:
                A[i+j*N_phi, i-1+j*N_phi] = 1/h_r**2 + 1/(2*(i+1)*h_r**2)  # i-1,j
            if j > 0:
                A[i+j*N_phi, i+(j-1)*N_phi] = 1/(h_r**2*h_phi**2*(i+1)**2)  # i,j-1
            if i < N_phi-1:
                A[i+j*N_phi, i+1+j*N_phi] = 1/h_r**2 + 1/(2*(i+1)*h_r**2)  # i+1,j
            if j < N_phi-1:
                A[i+j*N_phi, i+(j+1)*N_phi] = 1/(h_r**2*h_phi**2*(i+1)**2)  # i,j+1

    # problem last vrednosti
    t1=timer()
    w, X = eigh(A,B) #sym -> hitreje!
    t_eig=timer()-t1
    print('t eig posplošen = ', round(t_eig,3), ' s')
    t1=timer()
    ws, Xs = eigsh(A, k=24, M=B, which='LM', sigma=0) #sparse
    t_eigsh=timer()-t1
    print('t eigsh posplošen = ', round(t_eigsh,3), ' s')
    return (w,ws), (X,Xs), (t_eig,t_eigsh), (x,y)
#######################################################################################
# snovne lastnosti - HOMOGENA
gamma=1 # RP sila vpetja, N/m
d=1 #debelina, m
rho=1 #kg/m^3
ch=np.sqrt(gamma/rho/d) #hitrost valovanja, m/s
print('c = ',ch,' m/s')

N_r=N_phi=1+ 50 #liho!
(wh,wsh,wp), (Xh,Xsh,Xp), (t_h,t_sh,t_p), (x,y) = homog(N_r)

# homogena
Nc,Nr=6,4
l=2
fig=plt.figure(figsize=(Nc*l*2,Nr*l+l))
gs=fig.add_gridspec(Nr,Nc, hspace=0.11, wspace=0.11)
ax=gs.subplots(sharex=True, sharey=True)
vmax,vmin=np.max(np.abs(Xh[:,:Nr*Nc])), np.min(np.abs(Xh[:,:Nr*Nc]))
# vmax,vmin=np.max((Xh[:,:Nr*Nc])), np.min((Xh[:,:Nr*Nc])) #oscilatorno!!!
i=j=0
for k in range(Nr*Nc):
    ax[j,i].contourf(x,y,np.abs(Xh[:,k]).reshape(N_r,N_phi), cmap='turbo', vmax=vmax, vmin=vmin)
    # ax[j,i].contourf(x,y,(Xh[:,k]).reshape(N_r,N_phi), cmap='turbo', vmax=vmax, vmin=vmin)
    # ax[j,i].scatter(x,y,c=np.abs(Xh[:,k]), cmap='turbo', vmax=vmax, vmin=vmin, s=5)
    ax[j,i].grid(False)
    ax[j,i].set_xticks([])
    ax[j,i].set_yticks([])
    ax[j, i].set_frame_on(False)
    ax[j, i].set_aspect('equal')
    ax[j,i].set_title(f'$\lambda_{{{k+1}}}$', fontsize=12)
    # ax[j,i].set_title(f'$\lambda_{{{k+1}}} = {(np.sqrt(-wh[k])*ch/2/np.pi):<.2f} Hz$', fontsize=12)
    i+=1
    if i==Nc: i=0; j+=1
cbar_ax = fig.add_axes([0.86, 0, 0.07, 1])
cbar_ax.axis('off')
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
sm = cm.ScalarMappable(cmap='turbo', norm=norm)
sm.set_array([])
cbar=fig.colorbar(sm, ax=cbar_ax)
plt.savefig('./slike/homog')
plt.show()

# homogena sparse
fig=plt.figure(figsize=(Nc*l*2,Nr*l+l))
gs=fig.add_gridspec(Nr,Nc, hspace=0.11, wspace=0.11)
ax=gs.subplots(sharex=True, sharey=True)
vmax,vmin=np.max(np.abs(Xsh[:,:Nr*Nc])), np.min(np.abs(Xsh[:,:Nr*Nc]))
i=j=0
for k in range(Nr*Nc):
    ax[j,i].contourf(x,y,np.abs(Xsh[:,k]).reshape(N_r,N_phi), cmap='turbo', vmax=vmax, vmin=vmin)
    ax[j,i].grid(False)
    ax[j,i].set_xticks([])
    ax[j,i].set_yticks([])
    ax[j, i].set_frame_on(False)
    ax[j, i].set_aspect('equal')
    ax[j,i].set_title(f'$\lambda_{{{k+1}}}$', fontsize=12)
    # ax[j,i].set_title(f'$\lambda_{{{k+1}}} = {(np.sqrt(-wh[k])*ch/2/np.pi):<.2f} Hz$', fontsize=12)
    i+=1
    if i==Nc: i=0; j+=1
cbar_ax = fig.add_axes([0.86, 0, 0.07, 1])
cbar_ax.axis('off')
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
sm = cm.ScalarMappable(cmap='turbo', norm=norm)
sm.set_array([])
cbar=fig.colorbar(sm, ax=cbar_ax)
plt.savefig('./slike/homog_sparse')
plt.show()

#potenčna
plt.figure(figsize=(6,2.5))
plt.contourf(x,y,np.abs(Xp).reshape(N_r,N_phi), cmap='turbo')
plt.colorbar()
plt.gca().set_aspect('equal')
plt.gca().set_xticks([])
plt.gca().set_yticks([])
plt.gca().set_frame_on(False)
plt.grid(False)
plt.savefig('./slike/potenčna')
plt.show()

#########################################################################################
# NEHOMOGENA - radialno naraščanje/upadanje gostote
N_r=N_phi=1+ 50 #liho!
gamma=1 # RP sila vpetja, N/m
d=1 #debelina, m
rho_R=10 #zunaj kg/m^3
rho_0=1 #sredina kg/m^3
dr=(rho_R-rho_0)/(N_r-1)
rho=np.empty((N_r,N_phi))
for i in range(N_r):
    rho[i,:] = rho_0+i*dr
c=np.sqrt(gamma/rho/d) #hitrost valovanja, m/s

(w,ws), (X,Xs), (t_n,t_ns), (x,y) = nehomog(rho)

plt.figure(figsize=(6,2.5))
plt.contourf(x,y,rho.T, cmap='binary')
plt.colorbar(label='$\\rho$')
plt.gca().set_aspect('equal')
plt.grid(False)
plt.savefig('./slike/gostota')
plt.show()

#nehomogena
Nc,Nr=6,4
l=2
fig=plt.figure(figsize=(Nc*l*2,Nr*l+l))
gs=fig.add_gridspec(Nr,Nc, hspace=0.11, wspace=0.11)
ax=gs.subplots(sharex=True, sharey=True)
vmax,vmin=np.max(np.abs(X[:,:Nr*Nc])), np.min(np.abs(X[:,:Nr*Nc]))
i=j=0
for k in range(Nr*Nc):
    ax[j,i].contourf(x,y,np.abs(X[:,k]).reshape(N_r,N_phi), cmap='turbo', vmax=vmax, vmin=vmin)
    ax[j,i].grid(False)
    ax[j,i].set_xticks([])
    ax[j,i].set_yticks([])
    ax[j, i].set_frame_on(False)
    ax[j, i].set_aspect('equal')
    ax[j,i].set_title(f'$\lambda_{{{k+1}}}$', fontsize=12)
    # ax[j,i].set_title(f'$\lambda_{{{k+1}}} = {(np.sqrt(-w[k])*c/2/np.pi):<.2f} Hz$', fontsize=12)
    i+=1
    if i==Nc: i=0; j+=1
cbar_ax = fig.add_axes([0.86, 0, 0.07, 1])
cbar_ax.axis('off')
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
sm = cm.ScalarMappable(cmap='turbo', norm=norm)
sm.set_array([])
cbar=fig.colorbar(sm, ax=cbar_ax)
plt.savefig('./slike/nehomog')
plt.show()

#nehomogena sparse
fig=plt.figure(figsize=(Nc*l*2,Nr*l+l))
gs=fig.add_gridspec(Nr,Nc, hspace=0.11, wspace=0.11)
ax=gs.subplots(sharex=True, sharey=True)
vmax,vmin=np.max(np.abs(Xs[:,:Nr*Nc])), np.min(np.abs(Xs[:,:Nr*Nc]))
i=j=0
for k in range(Nr*Nc):
    ax[j,i].contourf(x,y,np.abs(Xs[:,k]).reshape(N_r,N_phi), cmap='turbo', vmax=vmax, vmin=vmin)
    ax[j,i].grid(False)
    ax[j,i].set_xticks([])
    ax[j,i].set_yticks([])
    ax[j, i].set_frame_on(False)
    ax[j, i].set_aspect('equal')
    ax[j,i].set_title(f'$\lambda_{{{k+1}}}$', fontsize=12)
    # ax[j,i].set_title(f'$\lambda_{{{k+1}}} = {(np.sqrt(-w[k])*c/2/np.pi):<.2f} Hz$', fontsize=12)
    i+=1
    if i==Nc: i=0; j+=1
cbar_ax = fig.add_axes([0.86, 0, 0.07, 1])
cbar_ax.axis('off')
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
sm = cm.ScalarMappable(cmap='turbo', norm=norm)
sm.set_array([])
cbar=fig.colorbar(sm, ax=cbar_ax)
plt.savefig('./slike/nehomog_sparse')
plt.show()

#lastne vrednosti
n=20
plt.plot(wh[:n], 'o', label='homogena')
plt.plot(wsh[:n], '1', label='homogena sparse')
plt.plot(w[:n], '+', label='nehomogena')
plt.plot(ws[:n], 'x', label='nehomogena sparse')
plt.axhline(wp, color='k', ls='--', label='potenčna')
plt.legend()
plt.xlabel('indeks')
plt.ylabel('$\lambda$ lastna vrednost')
plt.savefig('./slike/eigval')
plt.show()

#čas
plt.bar([1,2,3,4,5], [t_h,t_sh,t_p,t_n,t_ns], color='gray')
plt.xticks([1,2,3,4,5], ['homogena', 'homogena\nsparse', 'potenčna', 'nehomogena', 'nehomogena\nsparse'])
# plt.bar([1,2], [t_n,t_ns], color='gray')
# plt.xticks([1,2], ['nehomogena', 'nehomogena\nsparse'])
plt.ylabel('$t$ [s]')
plt.grid(False)
plt.savefig('./slike/čas')
plt.show()