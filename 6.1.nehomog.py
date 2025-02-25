#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 6. naloga: PDE lastne rešitve
# NEhomogena kvadratn opna

import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer
from numba import njit,jit
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from scipy.linalg import eigh,eig
from scipy.sparse.linalg import eigs,eigsh

plt.style.use('aljaz')
np.set_printoptions(linewidth=200, precision=3, suppress=True)

def g_homog(x,y,r):
    '''homogena opna'''
    return np.ones_like(r) #homogena opna
def g_kvad(x,y,r, ratio=2, sr=True):
    '''kvadratno območje odebelitve'''
    l = 0.5*a #stranica kvadrata
    if sr: origin = (a/2,b/2) #izhodišče
    else: origin = (3*a/4, 3*b/4) #desno zgoraj
    for i in range(r.shape[0]):
        for j in range(r.shape[1]):
            if abs(x[i]-origin[0]) < l/2 and abs(y[j]-origin[1]) < l/2:
                r[j,i] = ratio
            else: r[j,i] = 1
    return r #nehomogena - debelejša v kvadratu
def g_f(x,y,r, ratio=2):
    '''oblika f'''
    for i in range(r.shape[0]):
        for j in range(r.shape[1]):
            if a/4 < x[i] < a/2:
                r[j,i] = ratio
            elif y[j] > 3*b/4 and x[i] > a/4:
                r[j,i] = ratio
            elif b/4 < y[j] < b/2 and x[i] < 3*a/4:
                r[j,i] = ratio
            else: r[j,i] = 1
    return r
@njit(cache=True)
def potenčna(A, tol=1e-3, Nmax=10000):
    """potenčna metoda za iskanje največje lastne vrednosti matrike A"""
    x=np.random.rand(A.shape[0]) #poišči boljši približek!
    l=0
    l_old=1
    i=0
    while i<Nmax and abs(l-l_old)>tol:
        l_old=l
        y=A@x
        x=y/np.linalg.norm(y,2)
        l=x.T@A@x
        i+=1
    print('max lambda =', round(l,3), ', iteracij = ', i, ', tol = ', abs(l-l_old))
    return l, x
@njit(cache=True)
def inv_potenčna(A, sigma, tol=1e-9,Nmax=10000):
    """inverzna potenčna metoda za iskanje lastne vrednosti matrike A blizu sigma"""
    A_inv=np.linalg.inv(A-sigma*np.empty_like(A))
    x=np.random.rand(A.shape[0])
    l=0
    l_old=1
    i=0
    L=[]
    while i<Nmax and abs(l-l_old)>tol:
        l_old=l
        y=A_inv@x
        x=y/np.linalg.norm(y,2)
        l=x.T@A@x
        i+=1
        L.append(l)
    print('lambda =', round(l,3), ', iteracij = ', i, ', tol = ', abs(l-l_old))
    return l, x, L
# @jit(cache=True, nopython=False) #ni očitne izboljšave!
def run(N, rho, fast=False):
    # snovne lastnosti
    gamma=1 # RP sila vpetja, N/m
    d=1 #debelina, m

    # matrika A
    A = np.zeros((N**2,N**2))
    A1 = np.zeros_like(A)
    B = np.zeros_like(A) #desna stran -> lambda = omega**2
    for i in range(N):
        for j in range(N):
            B[i+j*N, i+j*N] = rho[i,j]*d/gamma
            A[i+j*N, i+j*N] = 4/h/h  # Main diagonal
            A1[i+j*N, i+j*N] = 4/h/h/rho[i,j]  # nesym problem! NI OK!!! TODO
            if i > 0:
                A[i+j*N, i-1+j*N] = -1/h/h  # Lower diagonal
                A1[i+j*N, i-1+j*N] = -1/h/h/rho[i-1,j]
            if j > 0:
                A[i+j*N, i+(j-1)*N] = -1/h/h  # Left diagonal
                A1[i+j*N, i+(j-1)*N] = -1/h/h/rho[i,j-1]
            if i < N-1:
                A[i+j*N, i+1+j*N] = -1/h/h  # Upper diagonal
                A1[i+j*N, i+1+j*N] = -1/h/h/rho[i+1,j]
            if j < N-1:
                A[i+j*N, i+(j+1)*N] = -1/h/h  # Right diagonal
                A1[i+j*N, i+(j+1)*N] = -1/h/h/rho[i,j+1]
    # plt.imshow(A1)
    # plt.show()

    # problem last vrednosti
    t1=timer()
    w, x = eigh(A,B) #scipy! posplošeni problem
    t_eig=timer()-t1
    print('t eig general = ', round(t_eig,3), ' s')
    if fast: return x, w, t_eig

    t1=timer()
    ws, xs = eigsh(A, k=9, M=B, which='SM') #sparse posplošeni
    t_eigs=timer()-t1
    ind=np.argsort(np.real(ws))
    ws=ws[ind]
    xs=xs[:,ind]
    print('t eig sparse = ', round(t_eigs,3), ' s')

    # return (x, xs), (w, ws), (t_eig, t_eigs)

    # t1=timer()
    # wp, xp = potenčna(A) #ne zanima nas največja lastna vred!, počasnejše
    # t_p=timer()-t1
    # print('t eig potenčna = ', round(t_p,3), ' s')

    t1=timer()
    wp, xp, L = inv_potenčna(A1, sigma=0) #najnižja lastna vrednost
    t_p=timer()-t1
    print('t eig sigma = ', round(t_p,3), ' s')
    # plt.plot(L,'kx-') #konvergenca
    # plt.xlabel('iteracije')
    # plt.ylabel('$\lambda_{min}$')
    # plt.show()
    
    t1=timer()
    wns, xns = np.linalg.eig(A1) #scipy! nesim, vendar klasični problem -> počasneje!
    # wns, xns = eig(A1) #scipy NI OK! OSCILACIJE! POČASNEJE!
    t_nesym=timer()-t1
    ind=np.argsort(np.real(wns))
    wns=wns[ind]
    xns=xns[:,ind]
    print('t eig nesym = ', round(t_nesym,3), ' s')

    t1=timer()
    wnss, xnss = eigs(A1, k=9, which='SM') #sparse nesym
    t_nesyms=timer()-t1
    ind=np.argsort(np.real(wnss))
    wnss=wnss[ind]
    xnss=xnss[:,ind]
    print('t eig sparse nesym = ', round(t_nesyms,3), ' s')

    return (x, xs, xns, xnss, xp), (w, ws, wns, wnss, wp), (t_eig, t_eigs, t_nesym, t_nesyms, t_p)
############################################################################################
# konvergenca f-opna

# N_list=np.append([3], np.linspace(10,70,7, dtype=int))
# # N_list=np.append([3], np.linspace(5,100,10*2, dtype=int))
# a=b=1 #širina=višina
# ratio = 100 #razmerje gostote

# # y=[1e-6,0.001,0.013,0.05,0.086,0.182,0.416,0.695,1.037,1.599,2.452,3.756,5.678,8.116,12.021,17.573,25.654,36.209,49.283,68.284,95.055]
# # plt.plot(N_list, y, 'kx-')
# # x=np.linspace(0,100,1000)
# # plt.plot(x, 1e-10*x**6, '-', c='gray')
# # plt.text(92, 92, '$N^6$', fontsize=12, color='k')
# # plt.xlabel('$N$')
# # plt.ylabel('$t$ [s]')
# # plt.minorticks_on()
# # plt.savefig('./slike/cas_N_new')
# # plt.show()

# čas=[]
# X=[]
# lamb=[]
# M=len(N_list)
# i=1
# for N in N_list:
#     print(f'{i}/{M}, N={N}', end='\t')
#     h=a/(N-1)
#     x = np.linspace(0,a,N)
#     y = np.linspace(0,b,N)
#     rho=np.zeros((len(x),len(y)))

#     rho = g_f(x,y,rho,ratio=ratio)

#     x,w,t = run(N, rho, fast=True)
#     čas.append(t)
#     X.append(x)
#     lamb.append(w)
#     i+=1

# plt.plot(N_list, čas, 'kx-')
# x=np.linspace(0,100,1000)
# plt.plot(x, 1e-10*x**6, '-', c='gray')
# plt.text(92, 92, '$N^6$', fontsize=12, color='k')
# plt.xlabel('$N$')
# plt.ylabel('$t$ [s]')
# plt.minorticks_on()
# # plt.savefig('./slike/cas_N')
# plt.show()

# plt.plot(N_list, [np.min(i) for i in lamb], 'kx-')
# plt.xlabel('$N$')
# plt.ylabel('$\lambda_{min}$')
# plt.minorticks_on()
# # plt.savefig('./slike/eigval_min_N')
# plt.show()

# plt.plot(N_list, [np.max(i) for i in X], 'kx-')
# plt.xlabel('$N$')
# plt.ylabel('$X_{max}$')
# plt.minorticks_on()
# # plt.savefig('./slike/Xmax_N')
# plt.show()

###################################################################################################
# # odvisnost od ratio
# r_list=np.append(1/np.logspace(0,2,10)[::-1],np.logspace(0,2,10)[1:])
# # r_list=np.append(1/np.logspace(0,2,20)[::-1],np.logspace(0,2,20)[1:])
# a=b=1 #širina=višina
# N = 10#50 #da je OK konvergenca!
# h=a/(N-1)
# x = np.linspace(0,a,N)
# y = np.linspace(0,b,N)
# rho=np.zeros((len(x),len(y)))

# čas={'f':[], 'sr':[], 'kot':[]}
# X={'f':[], 'sr':[], 'kot':[]}
# lamb={'f':[], 'sr':[], 'kot':[]}
# M=len(r_list)
# i=1
# for ratio in r_list:
#     print(f'{i}/{M}, ratio={ratio}', end='\t')

#     rho = g_f(x,y,rho,ratio=ratio)
#     xi,w,t = run(N, rho, fast=True)
#     čas['f'].append(t)
#     X['f'].append(xi)
#     lamb['f'].append(w)

#     rho = g_kvad(x,y,rho,ratio=ratio)
#     xi,w,t = run(N, rho, fast=True)
#     čas['sr'].append(t)
#     X['sr'].append(xi)
#     lamb['sr'].append(w)

#     rho = g_kvad(x,y,rho,ratio=ratio, sr=False)
#     xi,w,t = run(N, rho, fast=True)
#     čas['kot'].append(t)
#     X['kot'].append(xi)
#     lamb['kot'].append(w)
#     i+=1

# plt.plot(r_list, čas['f'], 'x-', label='f')
# plt.plot(r_list, čas['sr'], '.-', label='sredina')
# plt.plot(r_list, čas['kot'], '+-', label='kot')
# avg=np.average([čas['f']+čas['sr']+čas['kot']])
# plt.axhline(avg, ls='--', c='k')
# plt.xlabel('ratio')
# plt.ylabel('$t$ [s]')
# plt.minorticks_on()
# plt.tight_layout()
# plt.xscale('log')
# plt.legend(title='zgostitev opne')
# # plt.savefig('./slike/cas_ratio_')
# plt.show()

# plt.plot(r_list, [np.min(i) for i in lamb['f']], 'x-', label='f')
# plt.plot(r_list, [np.min(i) for i in lamb['sr']], '.-', label='sredina')
# plt.plot(r_list, [np.min(i) for i in lamb['kot']], '+-', label='vogal')
# plt.xlabel('ratio')
# plt.ylabel('$\lambda_{min}$')
# plt.xscale('log')
# plt.minorticks_on()
# plt.legend(title='zgostitev opne')
# # plt.savefig('./slike/eigval_min_ratio_')
# plt.show()

# plt.plot(r_list, [np.max(np.abs(i[:,0])) for i in X['f']], 'x-', label='f')
# plt.plot(r_list, [np.max(np.abs(i[:,0])) for i in X['sr']], '.-', label='sredina')
# plt.plot(r_list, [np.max(np.abs(i[:,0])) for i in X['kot']], '+-', label='vogal')
# plt.xlabel('ratio')
# plt.ylabel('1. lastna oblika $|X_{max}|$')
# plt.xscale('log')
# plt.legend(title='zgostitev opne')
# plt.minorticks_on()
# # plt.savefig('./slike/Xmax_ratio')
# plt.show()

##############################################################################################
# primerjava modalnih oblilk
def MAC(phi_l, phi_k):    
    """Modal Assurance Criterion (MAC)"""    
    MAC_mat = np.abs(np.einsum('ri,ik->rk',np.conj(phi_l).T,phi_k))**2 / \
            (np.einsum('ri,ir->r',np.conj(phi_l).T,phi_l)[:,np.newaxis] * \
             np.einsum('ri,ir->r',np.conj(phi_k).T,phi_k))
    return MAC_mat

# ratio=0.01
# a=b=1 #širina=višina
# N = 50 #da je OK konvergenca!
# h=a/(N-1)
# x = np.linspace(0,a,N)
# y = np.linspace(0,b,N)
# rho=np.zeros((len(x),len(y)))

# rho = g_f(x,y,rho,ratio=ratio)
# x1,w1,t1 = run(N, rho, fast=True)

# rho = g_kvad(x,y,rho,ratio=ratio)
# x2,w2,t2 = run(N, rho, fast=True)

# rho = g_kvad(x,y,rho,ratio=ratio, sr=False)
# x3,w3,t3 = run(N, rho, fast=True)

# n=10
# M12 = MAC(x1[:,:n], x2[:,:n])
# M13 = MAC(x1[:,:n], x3[:,:n])
# M23 = MAC(x2[:,:n], x3[:,:n])

# fig,ax=plt.subplots()
# plt.imshow(M12, origin='upper', aspect='equal', cmap='binary', vmin=0, vmax=1)
# plt.colorbar(ax=ax)
# plt.xlabel('$X$ f-opna')
# plt.ylabel('$X$ obtežba sredina')
# plt.grid(False)
# x=np.linspace(0,n-1,n,dtype=int)
# plt.xticks(x,x+1)
# plt.yticks(x,x+1)
# # plt.savefig(f'./slike/MAC_f_sr_r{ratio}.pdf')
# plt.show()

# fig,ax=plt.subplots()
# plt.imshow(M13, origin='upper', aspect='equal', cmap='binary', vmin=0, vmax=1)
# plt.colorbar(ax=ax)
# plt.xlabel('$X$ f-opna')
# plt.ylabel('$X$ obtežba vogal')
# plt.grid(False)
# plt.xticks(x,x+1)
# plt.yticks(x,x+1)
# # plt.savefig(f'./slike/MAC_f_vogal_r{ratio}.pdf')
# plt.show()

# fig,ax=plt.subplots()
# plt.imshow(M23, origin='upper', aspect='equal', cmap='binary', vmin=0, vmax=1)
# plt.colorbar(ax=ax)
# plt.xlabel('$X$ obtežba sredina')
# plt.ylabel('$X$ obtežba vogal')
# plt.grid(False)
# plt.xticks(x,x+1)
# plt.yticks(x,x+1)
# # plt.savefig(f'./slike/MAC_sr_vogal_r{ratio}.pdf')
# plt.show()

##################################################################################################
# prikaz
N=30 #št. računskih točk po osi
a=b=1 #širina=višina
ratio = 100 #razmerje gostote -> povzroča težave sparse solverju!?
h=a/(N-1)

x = np.linspace(0,a,N)
y = np.linspace(0,b,N)
rho=np.zeros((len(x),len(y)))

# rho = g_homog(x,y,rho)
# rho = g_kvad(x,y,rho,ratio=ratio, sr=False)
rho = g_f(x,y,rho,ratio=ratio)

# plt.imshow(rho, origin='lower', aspect='equal', cmap='binary')
# plt.xticks([])
# plt.yticks([])
# plt.savefig('./slike/opna')
# plt.show()

# c=np.sqrt(gamma/rho/d) #hitrost valovanja, m/s
# print('c_max = ',np.max(c), 'm/s , c_min = ',np.min(c),' m/s')

Nc,Nr=3,3
l=2
(x, xs, xns, xnss, xp), (w, ws, wns, wnss, wp), (t_eig, t_eigs, t_nesym, t_nesyms, t_p) = run(N, rho) #nesym NI OK, potenčna le homogena!
# (x, xs), (w, ws), (t_eig, t_eigs) = run(N, rho)

n=20
plt.plot(w[:n], 'o', label='posplošeno')
plt.plot(np.real(ws[:n]), 'x', label='sparse posplošeno')
plt.plot(np.real(wns[:n]), '+', label='ne-sym') #sortiranje!
plt.plot(np.real(wnss[:n]), '1', label='sparse ne-sym')
plt.axhline(wp, ls='--', c='red', label='inverzna potenčna\nmetoda ($\sigma=0$)') #le homog
plt.legend()
plt.xlabel('indeks')
plt.ylabel('$\lambda$ lastna vrednost')
plt.savefig('./slike/eigval')
plt.show()

# #inverzna potenčna metoda
# plt.imshow(xp.reshape(N,N), origin='lower', extent=(0,a,0,b), aspect='equal')
# plt.colorbar()
# plt.savefig('./slike/potencna_nehomog')
# plt.show()

# n=10 #OK korelacije!
# M12 = np.abs(MAC(x[:,:n], xs[:,:n]))
# M13 = np.abs(MAC(x[:,:n], xns[:,:n]))
# M23 = np.abs(MAC(x[:,:n], xnss[:,:n]))
# fig,ax=plt.subplots()
# plt.imshow(M12, origin='upper', aspect='equal', cmap='binary', vmin=0, vmax=1)
# plt.colorbar(ax=ax)
# plt.xlabel('$X$ posplošeno')
# plt.ylabel('$X$ sparse posplošeno')
# plt.grid(False)
# xx=np.linspace(0,n-1,n,dtype=int)
# plt.xticks(xx,xx+1)
# plt.yticks(xx,xx+1)
# # plt.savefig('./slike/')
# plt.show()

# fig,ax=plt.subplots()
# plt.imshow(M13, origin='upper', aspect='equal', cmap='binary', vmin=0, vmax=1)
# plt.colorbar(ax=ax)
# plt.xlabel('$X$ posplošeno')
# plt.ylabel('$X$ nesim')
# plt.grid(False)
# plt.xticks(xx,xx+1)
# plt.yticks(xx,xx+1)
# # plt.savefig('./slike/')
# plt.show()

# fig,ax=plt.subplots()
# plt.imshow(M23, origin='upper', aspect='equal', cmap='binary', vmin=0, vmax=1)
# plt.colorbar(ax=ax)
# plt.xlabel('$X$ posplošeno')
# plt.ylabel('$X$ sparse nesim')
# plt.grid(False)
# plt.xticks(xx,xx+1)
# plt.yticks(xx,xx+1)
# # plt.savefig('./slike/')
# plt.show()

#general
fig=plt.figure(figsize=(Nc*l,Nr*l+0.9))
gs=fig.add_gridspec(Nr,Nc, hspace=0.11, wspace=0.11)
ax=gs.subplots(sharex=True, sharey=True)
vmax,vmin=np.max((x[:,:Nr*Nc])), np.min((x[:,:Nr*Nc]))
i=j=0
for k in range(Nr*Nc):
    ax[j,i].imshow((x[:,k]).reshape(N,N), origin='lower', extent=(0,a,0,b), aspect='equal', vmax=vmax, vmin=vmin)
    ax[j,i].grid(False)
    ax[j,i].set_xticks([])
    ax[j,i].set_yticks([])
    ax[j,i].set_title(f'$\lambda_{{{k+1}}}$', fontsize=9)
    # ax[j,i].set_title(f'$\lambda_{{{k+1}}} = {(np.sqrt(-w[k])/2/np.pi):<.2f} Hz$', fontsize=9)
    i+=1
    if i==Nc: i=0; j+=1
cbar_ax = fig.add_axes([0.86, 0, 0.07, 1])
cbar_ax.axis('off')
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
sm = cm.ScalarMappable(cmap='turbo', norm=norm)
sm.set_array([])
cbar=fig.colorbar(sm, ax=cbar_ax)
cbar.ax.tick_params(labelsize=9)
plt.savefig('./slike/lastne')
plt.show()

#sparse general
fig=plt.figure(figsize=(Nc*l,Nr*l+0.9))
gs=fig.add_gridspec(Nr,Nc, hspace=0.11, wspace=0.11)
ax=gs.subplots(sharex=True, sharey=True)
vmax,vmin=np.max(np.abs(xs[:,:Nr*Nc])), np.min(np.abs(xs[:,:Nr*Nc]))
vmax,vmin=np.max(np.real(xs[:,:Nr*Nc])), np.min(np.real(xs[:,:Nr*Nc]))
i=j=0
for k in range(Nr*Nc):
    # ax[j,i].imshow(np.abs(xs[:,k]).reshape(N,N), origin='lower', extent=(0,a,0,b), aspect='equal', vmax=vmax, vmin=vmin)
    ax[j,i].imshow(np.real(xs[:,k]).reshape(N,N), origin='lower', extent=(0,a,0,b), aspect='equal', vmax=vmax, vmin=vmin)
    ax[j,i].grid(False)
    ax[j,i].set_xticks([])
    ax[j,i].set_yticks([])
    ax[j,i].set_title(f'$\lambda_{{{k+1}}}$', fontsize=9)
    # ax[j,i].set_title(f'$\lambda_{{{k+1}}} = {(np.sqrt(-ws[k])/2/np.pi):<.2f} Hz$', fontsize=9)
    i+=1
    if i==Nc: i=0; j+=1
cbar_ax = fig.add_axes([0.86, 0, 0.07, 1])
cbar_ax.axis('off')
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
sm = cm.ScalarMappable(cmap='turbo', norm=norm)
sm.set_array([])
cbar=fig.colorbar(sm, ax=cbar_ax)
cbar.ax.tick_params(labelsize=9)
plt.savefig('./slike/sparse')
plt.show()

#nesym
fig=plt.figure(figsize=(Nc*l,Nr*l+0.9))
gs=fig.add_gridspec(Nr,Nc, hspace=0.11, wspace=0.11)
ax=gs.subplots(sharex=True, sharey=True)
vmax,vmin=np.max(np.real(xns[:,:Nr*Nc])), np.min(np.real(xns[:,:Nr*Nc]))
i=j=0
for k in range(Nr*Nc):
    ax[j,i].imshow(np.real(xns[:,k]).reshape(N,N), origin='lower', extent=(0,a,0,b), aspect='equal', vmax=vmax, vmin=vmin)
    ax[j,i].grid(False)
    ax[j,i].set_xticks([])
    ax[j,i].set_yticks([])
    ax[j,i].set_title(f'$\lambda_{{{k+1}}}$', fontsize=9)
    # ax[j,i].set_title(f'$\lambda_{{{k+1}}} = {(np.sqrt(-ws[k])/2/np.pi):<.2f} Hz$', fontsize=9)
    i+=1
    if i==Nc: i=0; j+=1
cbar_ax = fig.add_axes([0.86, 0, 0.07, 1])
cbar_ax.axis('off')
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
sm = cm.ScalarMappable(cmap='turbo', norm=norm)
sm.set_array([])
cbar=fig.colorbar(sm, ax=cbar_ax)
cbar.ax.tick_params(labelsize=9)
plt.savefig('./slike/nesym')
plt.show()

#sparse nesym
fig=plt.figure(figsize=(Nc*l,Nr*l+0.9))
gs=fig.add_gridspec(Nr,Nc, hspace=0.11, wspace=0.11)
ax=gs.subplots(sharex=True, sharey=True)
vmax,vmin=np.max(np.real(xnss[:,:Nr*Nc])), np.min(np.real(xnss[:,:Nr*Nc]))
i=j=0
for k in range(Nr*Nc):
    ax[j,i].imshow(np.real(xnss[:,k]).reshape(N,N), origin='lower', extent=(0,a,0,b), aspect='equal', vmax=vmax, vmin=vmin)
    ax[j,i].grid(False)
    ax[j,i].set_xticks([])
    ax[j,i].set_yticks([])
    ax[j,i].set_title(f'$\lambda_{{{k+1}}}$', fontsize=9)
    # ax[j,i].set_title(f'$\lambda_{{{k+1}}} = {(np.sqrt(-ws[k])/2/np.pi):<.2f} Hz$', fontsize=9)
    i+=1
    if i==Nc: i=0; j+=1
cbar_ax = fig.add_axes([0.86, 0, 0.07, 1])
cbar_ax.axis('off')
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
sm = cm.ScalarMappable(cmap='turbo', norm=norm)
sm.set_array([])
cbar=fig.colorbar(sm, ax=cbar_ax)
cbar.ax.tick_params(labelsize=9)
plt.savefig('./slike/sparse_nesym')
plt.show()