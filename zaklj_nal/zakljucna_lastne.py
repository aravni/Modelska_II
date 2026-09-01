#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Zaključna naloga: upogib nosilca, lastne oblike

import numpy as np 
import matplotlib.pyplot as plt
from scipy.linalg import eig, eigh
from scipy.sparse.linalg import eigs, eigsh
import matplotlib.animation as animation

plt.style.use('aljaz')

N = 100 # št. vozlišč
Nps = 2 #izbira KE (2-nosilec, 3-okvir)

# snovne lastnosti
ρ = 1 # kg/mm^3
E = 100 # MPa
ν = 0.3 # Poisson

# gemeotrija nosilca
L = 100 # mm
alpha = np.deg2rad(0) # kot zasuka
RP = [1,1] # vpetje (0-prosto, 1-fiksno)

# geometrijske lastnosti prereza
# 1. pravokotnik: a-širina, h-višina
a, h = 10, 10 #mm
A = a*h #mm^2
I = 1/12 * a*h**3 #mm^4
# 2. krog: r-polmer
# r = 1 #mm
# A = np.pi*r**2 #mm^2
# I = 1/64 * np.pi*r**4 #mm^4

###################################################################################

def T2(α):
    """transformacijska matrika za 2 PS"""    
    T = np.array([[ np.cos(α), np.sin(α),         0,         0],
                  [-np.sin(α), np.cos(α),         0,         0],
                  [         0,         0, np.cos(α), np.sin(α)],
                  [         0,         0,-np.sin(α), np.cos(α)]])
    return T
def T3(α):
    """transformacijska matrika za 3 PS"""    
    T = np.array([[ np.cos(α), np.sin(α), 0,         0,         0, 0],
                  [-np.sin(α), np.cos(α), 0,         0,         0, 0],
                  [         0,         0, 1,         0,         0, 0],
                  [         0,         0, 0, np.cos(α), np.sin(α), 0],
                  [         0,         0, 0,-np.sin(α), np.cos(α), 0],
                  [         0,         0, 0,         0,         0, 1]])
    return T
def M_e_f(A, L, ρ, α):
    """masna matrika za KE okvir/frame"""
    M_lok = ρ*A*L*np.array([[1/3,         0,         0, 1/6,         0,         0],
                            [  0,     13/35,  11*L/210,   0,      9/70, -13*L/420],
                            [  0,  11*L/210,  L**2/105,   0,  13*L/420, -L**2/140],
                            [1/6,         0,         0, 1/3,         0,         0],
                            [  0,      9/70,  13*L/420,   0,     13/35, -11*L/210],
                            [  0, -13*L/420, -L**2/140,   0, -11*L/210,  L**2/105]])
    return T3(α).T @ M_lok @ T3(α)
def K_e_f(A, E, I, L, α):
    """togostna matrika za KE okvir/frame"""
    K_lok = E/L*np.array([[ A,          0,      0, -A,          0,      0],
                          [ 0,  12*I/L**2,  6*I/L,  0, -12*I/L**2,  6*I/L],
                          [ 0,      6*I/L,    4*I,  0,     -6*I/L,    2*I],
                          [-A,          0,      0,  A,          0,      0],
                          [ 0, -12*I/L**2, -6*I/L,  0,  12*I/L**2, -6*I/L],
                          [ 0,      6*I/L,    2*I,  0,     -6*I/L,    4*I]])
    return T3(α).T @ K_lok @ T3(α)
def M_e_n(A, Le, ρ, α):
    """masna matrika za KE Euler-Bernoulli nosilec"""
    M_lok = ρ*A*Le/420*np.array([[156,     22*Le,    54,    -13*Le],
                                 [22*Le,   4*Le**2,  13*Le, -3*Le**2],
                                 [54,      13*Le,    156,   -22*Le],
                                 [-13*Le, -3*Le**2, -22*Le,  4*Le**2]])
    return T2(α).T @ M_lok @ T2(α)
def K_e_n(E, I, Le, α):
    """togostna matrika za KE Euler-Bernoulli nosilec"""
    K_lok = E*I/Le**3*np.array([[12,   6*Le,   -12,   6*Le],
                                [6*Le, 4*Le**2,-6*Le, 2*Le**2],
                                [-12, -6*Le,    12,  -6*Le],
                                [6*Le, 2*Le**2,-6*Le, 4*Le**2]])
    return T2(α).T @ K_lok @ T2(α)

###################################################################################

# vozlišča in elementi
vozlisca0 = np.transpose([np.linspace(0,L,N), np.zeros(N)])
vozlisca = np.empty_like(vozlisca0)
vozlisca[:,0] = np.cos(alpha)*vozlisca0[:,0] - np.sin(alpha)*vozlisca0[:,1]
vozlisca[:,1] = np.sin(alpha)*vozlisca0[:,0] + np.cos(alpha)*vozlisca0[:,1]
ind = np.linspace(0,N-1,N, dtype=int)
elementi = np.transpose([ind[:-1], ind[1:]]).astype(int)

# sestav globalne masne in togostne matrike
M_glob = np.zeros((N*Nps,N*Nps))
K_glob = np.zeros((N*Nps,N*Nps))

def d(a,b): return ((a[0]-b[0])**2 + (a[1]-b[1])**2)**0.5
for element in elementi:
    Le = d(vozlisca[element[0]], vozlisca[element[1]]) # dolžina elementa
    αe = np.arctan2(np.diff(vozlisca[element,1]) , np.diff(vozlisca[element,0]))[0] # kot zasuka
    # priševanje prispevkov posameznih elementov
    if Nps == 3: #frame/okvir
        ind = np.array([[Nps*element[0],Nps*element[0]+1,Nps*element[0]+2,Nps*element[1],Nps*element[1]+1,Nps*element[1]+2]])
        M_glob[ind.T,ind] += M_e_f(A, Le, ρ, αe)
        K_glob[ind.T,ind] += K_e_f(A, E, I, Le, αe)
    elif Nps == 2: #nosilec
        ind = np.array([[Nps*element[0],Nps*element[0]+1,Nps*element[1],Nps*element[1]+1]])
        M_glob[ind.T,ind] += M_e_n(A, Le, ρ, αe)
        K_glob[ind.T,ind] += K_e_n(E, I, Le, αe)
    
# robni pogoji
ind = np.array([np.arange(N*Nps)],dtype=int)
if RP[0] == 1: #prvo vozlisce je fiksno
    ind = ind[:,Nps:]
if RP[1] == 1: #zadnje vozlisce je fiksno
    ind = ind[:,:-Nps]
M_glob = M_glob[ind.T,ind]
K_glob = K_glob[ind.T,ind]

# lastne vrednosti in lastni vektorji
eig_val, eig_vec = eig(K_glob, M_glob)

# urejanje po velikosti
eig_ind = np.argsort(np.abs(eig_val))
eig_val = eig_val[eig_ind]
eig_vec = eig_vec[:,eig_ind]
eig_freq = np.sqrt(np.abs(eig_val)) / 2 / np.pi # lastne frekvence [Hz]

# modalne oblike razširimo še na prostostne stopnje kjer ni pomikov
eig_vec_ = np.zeros((N*Nps,eig_vec.shape[0]))
eig_vec_[ind,:] = eig_vec

if Nps == 3: # izbrišemo zasuke - vsaka tretja PS (eig_vec_[::3,:])
    eig_vec_ = eig_vec_[(np.arange(eig_vec_.shape[0])+1) % 3 != 0, :]

########################################################################################

# izris lastne oblike
N_max = 7 # št. lastnih oblik za izris
m = 0 # zaporedni indeks lastne oblike
scale = 1 # prečno nihanje
scale_x = 0.2 # vzdolžno nihanje

def animate(i, smer, line):
    if smer == 'x':
        # return ax.plot(amp[i]*modeshape[:,0], modeshape[:,1],'k-o')
        line.set_data(amp[i]*modeshape[:,0], modeshape[:,1])
    else:
        # return ax.plot(modeshape[:,0], amp[i]*modeshape[:,1],'k-o')
        line.set_data(modeshape[:,0], amp[i]*modeshape[:,1])
    return line,
        
for m in range(N_max):
    modeshape = np.copy(vozlisca).ravel()
    modeshape += eig_vec_[:,m]
    modeshape = modeshape.reshape(N,2)

    print(f'mode {m+1}, f = {eig_freq[m]:.5g} Hz')

    # max amplituda
    plt.clf()
    # nedeformirano stanje
    [plt.plot(vozlisca[[i,j],0],vozlisca[[i,j],1],'-', c='gray') for i,j in elementi]
    plt.plot(vozlisca[:,0],vozlisca[:,1],'o', c='gray')
    # deformirano stanje
    [plt.plot(modeshape[[i,j],0], modeshape[[i,j],1],'-', c='k') for i,j in elementi]
    plt.plot(modeshape[:,0], modeshape[:,1],'o', c='k')
    [plt.plot(modeshape[[i,j],0], -modeshape[[i,j],1],'-', c='k') for i,j in elementi]
    plt.plot(modeshape[:,0], -modeshape[:,1],'o', c='k')
    # plt.axis('equal')
    plt.xlabel('$x$ [mm]')
    plt.ylabel('$y$ [mm]')
    plt.minorticks_on()
    plt.tight_layout()
    # plt.savefig(f'./slike/dinamika/prostoD_mode_{m+1}.pdf')
    plt.show()

    # animacija
    plt.clf()
    fig,ax = plt.subplots()
    [ax.plot(vozlisca[[i,j],0],vozlisca[[i,j],1],'-', c='gray') for i,j in elementi]
    ax.plot(vozlisca[:,0],vozlisca[:,1],'o', c='gray')
    ax.set_xlabel('$x$ [mm]')
    ax.set_ylabel('$y$ [mm]')
    ax.set_title(f'mode ${m+1}$, $f={eig_freq[m]:.3g}$ Hz')
    plt.minorticks_on()
    plt.tight_layout()

    if np.max(modeshape[:,1]) != 0:
        smer, s = 'y', scale
        ymax=max([np.max(modeshape[:,1]),np.max(vozlisca[:,1])]) 
        ax.set_xlim(0, max([np.max(modeshape[:,0]), np.max(vozlisca[:,0])])*s*1.1)
        ax.set_ylim(-ymax*s*1.1, ymax*s*1.1)

        amp = scale*np.linspace(-1,1,100)
        amp = np.block([amp,amp[::-1]])

        line, = ax.plot(modeshape[:,0], amp[0]*modeshape[:,1],'k-o')
    else:
        smer, s = 'x', scale_x #ni nihanja v y smeri
        ax.set_xlim(0, max([np.max(modeshape[:,0]), np.max(vozlisca[:,0])])*(1+s)*1.1)

        amp = np.linspace(1-s,1+s,100)
        amp = np.block([amp,amp[::-1]])

        line, = ax.plot(amp[0]*modeshape[:,0], modeshape[:,1],'k-o')

    ani = animation.FuncAnimation(
        fig, animate, interval=5, blit=True, repeat=True, frames=len(amp), fargs=(smer, line))
    # ani.save(f'./slike/dinamika/mode_{m+1}.gif')
    # ani.save(f'./slike/dinamika/prosto_mode_{m+1}.mp4')
    plt.show()
    plt.close(fig)