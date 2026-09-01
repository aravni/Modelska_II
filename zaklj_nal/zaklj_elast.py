#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Zaključna naloga: upogib nosilca, statika
# linearno elastična reološka zveza, Euler-Bernoulli in Timošenkov nosilec
# vir: https://vulcanhammer.net/wp-content/uploads/2017/01/finite-elements-in-plasticity-theory-and-practice.pdf

import numpy as np 
import matplotlib.pyplot as plt
from scipy import linalg
from scipy.sparse.linalg import spsolve
from timeit import default_timer as timer
import matplotlib
from alignAxes import *

plt.style.use('aljaz')

def K_e_eb(E, I, Le):
    """togostna matrika za KE Euler-Bernoulli nosilec"""
    K_lok = E*I/Le**3*np.array([[12,   6*Le,   -12,   6*Le],
                                [6*Le, 4*Le**2,-6*Le, 2*Le**2],
                                [-12, -6*Le,    12,  -6*Le],
                                [6*Le, 2*Le**2,-6*Le, 4*Le**2]])
    return K_lok
def K_e_tg2(E, I, Le):
    """togostna matrika za KE Timošenkov nosilec Gaussova integracija v 2 točkah
    - pretoge rešitve po [https://vulcanhammer.net/wp-content/uploads/2017/01/finite-elements-in-plasticity-theory-and-practice.pdf]"""
    M1 = E*I/Le*np.array([[0, 0, 0, 0],
                          [0, 1, 0,-1],
                          [0, 0, 0, 0],
                          [0,-1, 0, 1]])
    M2 = G*As/Le*np.array([[1,     Le/2,   -1,     Le/2],
                           [Le/2,  Le**2/3,-Le/2,  Le**2/6],
                           [-1,   -Le/2,    1,    -Le/2],
                           [Le/2,  Le**2/6,-Le/2,  Le**2/3]])
    return M1+M2
def K_e_tg1(E, I, Le):
    """togostna matrika za KE Timošenkov nosilec Gaussova integracija v 1 točki"""
    M1 = E*I/Le*np.array([[0, 0, 0, 0],
                          [0, 1, 0,-1],
                          [0, 0, 0, 0],
                          [0,-1, 0, 1]])
    M2 = G*As/Le*np.array([[1,     Le/2,   -1,     Le/2],
                           [Le/2,  Le**2/4,-Le/2,  Le**2/4],
                           [-1,   -Le/2,    1,    -Le/2],
                           [Le/2,  Le**2/4,-Le/2,  Le**2/4]])
    return M1+M2
def K_e_tm(E, I, Le):
    """togostna matrika za KE Timošenkov nosilec modificirana oblika
    - kubična polinomska aproksimacija povesa in zasuka -> OK!"""
    C = (12*E*I)/(As*G*Le**2)
    K_lok = E*I/Le**3/(1+C)*np.array([[12,   6*Le,       -12,   6*Le],
                                      [6*Le, (4+C)*Le**2,-6*Le, (2-C)*Le**2],
                                      [-12, -6*Le,        12,  -6*Le],
                                      [6*Le, (2-C)*Le**2,-6*Le, (4+C)*Le**2]])
    return K_lok

N = 1000 # št. vozlišč
KE = 'EB' #izbira KE (EB=nosilec Euler-Bernoulli, TG1,TG2,TM-oblike za Timošenko)
# KE = 'EB' #izbira KE (EB=nosilec Euler-Bernoulli, TG1,TG2,TM-oblike za Timošenko)

# snovne lastnosti
R_p = 1 # meja plastičnosti, MPa
ρ = 1 # kg/mm^3
E = 1 # MPa
ν = 0.3 # Poisson
G = E/(2*(1+ν)) # strižni modul v elast. območju, MPa

# gemeotrija nosilca
L = 100 # mm
RP = [2,0] # vpetje (0-prosto, 1-fiksno, 2-konzolno) [levo, desno]
RP_F = [0,0,1] # točkovna sila F0 (0-neobremenjeno, 1-obremenjeno) [levo, sredina, desno]
RP_M = [0,0,0] # upogibni moment M0 (0-neobremenjeno, 1-obremenjeno) [levo, sredina, desno]
F0 = 100 # točkovna sila - navzdol, N
M0 = 0e3 # upogibni moment - smer ure, N mm
q0 = 0 # linijska obremenitev, N/mm

# geometrijske lastnosti prereza
# # 1. pravokotnik: a-širina, h-višina
a, h = 100, 100 #mm
κ = 10*(1+ν)/(12+11*ν) #5/6 #koef. oblike prereza
A = a*h #mm^2
I = 1/12 * a*h**3 #mm^4
# 2. krog: r-polmer
# r = 100 #mm
# κ = 6*(1+ν)/(7+6*ν) #9/10 #koef. oblike prereza
# A = np.pi*r**2 #mm^2
# I = 1/64 * np.pi*r**4 #mm^4

hx = L/(N-1)
As = κ*A #vpliv oblike prereza na velikost max strižnih nap., upoštevajoč konst. porazd. stižnih nap.

# reševanje sistema enačb za pomike
#############################################################################

# vozlišča in elementi
vozlisca = np.transpose([np.linspace(0,L,N), np.zeros(N)])
ind = np.linspace(0,N-1,N, dtype=int)
elementi = np.transpose([ind[:-1], ind[1:]]).astype(int)
N_e = len(elementi) #=N-1

# sestav globalne masne in togostne matrike
Nps = 2
M_glob = np.zeros((N*Nps,N*Nps))
K_glob = np.zeros((N*Nps,N*Nps))
F_glob = np.zeros(N*Nps) # desna stran - vozliščne obremenitve

def d(a,b): return ((a[0]-b[0])**2 + (a[1]-b[1])**2)**0.5
for element in elementi:
    Le = d(vozlisca[element[0]], vozlisca[element[1]]) # dolžina elementa (Le = hx)
    # priševanje prispevkov posameznih elementov
    ind = np.array([[Nps*element[0],Nps*element[0]+1,Nps*element[1],Nps*element[1]+1]])
    if KE == 'EB': K_glob[ind.T,ind] += K_e_eb(E, I, Le)
    elif KE == 'TG1': K_glob[ind.T,ind] += K_e_tg1(E,I,Le)
    elif KE == 'TG2': K_glob[ind.T,ind] += K_e_tg2(E,I,Le)
    elif KE == 'TM': K_glob[ind.T,ind] += K_e_tm(E,I,Le)
    else: raise Exception('Napačen tip KE!')
    F_glob[ind[0,[0,2]]] += -q0*Le/2
# plt.imshow(np.abs(K_glob), cmap='binary')
# plt.colorbar()
# plt.show()

# robni pogoji
if RP[0] == 1 or RP[0] == 2: #prvo vozlisce je fiksno
    i = 0
    K_glob[i,:] = 0
    K_glob[i,i] = 1
    F_glob[i] = 0
if RP[0] == 2: #prvo vozlisce konzola
    i = 1
    K_glob[i,:] = 0
    K_glob[i,i] = 1
    F_glob[i] = 0
if RP[1] == 1 or RP[1] == 2: #zadnje vozlisce je fiksno
    i = -2
    K_glob[i,:] = 0
    K_glob[i,i] = 1
    F_glob[i] = 0
if RP[1] == 2: #zadnje vozlisce konzola
    i = -1
    K_glob[i,:] = 0
    K_glob[i,i] = 1
    F_glob[i] = 0

# obremenitve
if RP_F[0] == 1: #prvo vozlisce sila F0
    F_glob[0] += -F0
if RP_F[1] == 1: #središčno vozlisce sila F0
    F_glob[(N*Nps)//2] += -F0
if RP_F[2] == 1: #zadnje vozlisce sila F0
    F_glob[-2] += -F0

if RP_M[0] == 1: #prvo vozlisce moment 
    F_glob[1] += -M0
if RP_M[1] == 1: #središčno vozlisce moment
    F_glob[(N*Nps)//2+1] += -M0
if RP_M[2] == 1: #zadnje vozlisce moment
    F_glob[-1] += -M0

# reševanje sistema enačb
t1=timer()
u1 = linalg.solve(K_glob, F_glob)
t2=timer()
u=u2 = spsolve(K_glob, F_glob)
t3=timer()
u = np.linalg.solve(K_glob, F_glob)
t4=timer()
U, Φ = u.reshape(vozlisca.shape).T
print(f'norme rešitev: scipy {np.linalg.norm(u1):.8g}, sparse {np.linalg.norm(u2):.8g}, numpy {np.linalg.norm(u):.8g}')
print(f'časi reševanja: scipy {(t2-t1):.3g} s, sparse {(t3-t2):.3g} s, numpy {(t4-t3):.3g} s')
print(f'maksimalni poves: {np.max(np.abs(U)):.8g} mm')

# analitična rešitev
x=np.linspace(0,L,N)
if KE == 'EB':
    u_max_ex = -F0*L**3/3/E/I #EB
    u_ex = -F0*x**2/6/E/I*(3*L-x) #Euler-Bernoulli
    φ_ex = -F0/E/I*(L*x-x**2/2) #=du/dx EB?
else:
    u_ex = -F0/As/G*x - F0/E/I*(-x**3/6+x**2*L/2) #Timošenko
    u_max= -F0*L**3/3/E/I -F0/As/G
    φ_ex = -F0/E/I*(L*x-x**2/2) - F0/As/G #=du/dx?
print(f'analitično max poves: Euler-Bernoulli {F0*L**3/3/E/I:.5g} mm, Timošenko {F0/As/G*L + F0/E/I*(-L**3/6+L**2*L/2):.5g} mm')
print(f'error u    : {np.linalg.norm(u_ex-U)/np.linalg.norm(u_ex):.3g}, {np.linalg.norm(u_ex-U):.3g}')
print(f'error phi  : {np.linalg.norm(φ_ex-Φ)/np.linalg.norm(φ_ex):.3g}, {np.linalg.norm(φ_ex-Φ):.3g}')
print(f'error u_max: {abs(np.max(np.abs(u_ex))-np.max(np.abs(U))):.3g}, {abs(np.max(np.abs(φ_ex))-np.max(np.abs(Φ))):.3g}')
input("Press Enter to continue...")

# izris
# nedeformirano stanje
scale=1
[plt.plot(vozlisca[[i,j],0], vozlisca[[i,j],1],'-', c='gray') for i,j in elementi]
plt.plot(vozlisca[:,0], vozlisca[:,1],'o', c='gray')
# deformirano stanje
[plt.plot(vozlisca[[i,j],0], U[[i,j]]*scale,'-', c='k') for i,j in elementi]
plt.plot(vozlisca[:,0], U*scale,'o', c='k')
plt.plot(x, u_ex*scale, 'r-')
plt.axis('equal')
plt.xlabel('$x$ [mm]')
plt.ylabel('$y$ [mm]')
plt.minorticks_on()
plt.tight_layout()
plt.show()

# plt.plot(vozlisca[:,0], Φ, 'ko-', label='rešitev')
# plt.plot(x, φ_ex, 'r--', label='analitično')
# plt.show()
# φ_ex = -F0/E/I*(L*vozlisca[:,0]-vozlisca[:,0]**2/2) - F0/As/G #=du/dx?
# φ_exb = -F0/E/I*(L*vozlisca[:,0]-vozlisca[:,0]**2/2) #=du/dx?
# phi =  np.gradient(U, hx, edge_order=2) #to ne velja za Timošenko!!! -> a enako analitično?
# print('phi',np.sum(np.abs(Φ-φ_ex)), np.sum(np.abs(Φ-φ_exb)), np.sum(np.abs(phi-φ_ex)), np.sum(np.abs(φ_exb-phi)))
# plt.plot(vozlisca[:,0], Φ-φ_ex, '-', label='timošenko')
# plt.plot(vozlisca[:,0], Φ-φ_exb, '-', label='EB')
# plt.show()

# post process
###########################################################################################
#TODO nekaj ni OK!! - premajhne vrednosti za EB
def epsilon_x_mkr(z):
    if KE == 'EB': #Euler-Bernoulli
        ε_xx = np.empty(N) #deformacija je drugi odvod pomikov, 2. red natančnosti
        ε_xx[0] = z * (2*U[0] - 5*U[1] + 4*U[2] - U[3]) / abs(vozlisca[2,0] - vozlisca[1,0])*abs(vozlisca[1,0] - vozlisca[0,0])
        ε_xx[-1] = z * (2*U[-1] - 5*U[-2] + 4*U[-3] - U[-4]) / abs(vozlisca[-3,0] - vozlisca[-2,0])*abs(vozlisca[-2,0] - vozlisca[-1,0])
        for i in range(1,N-1):
            ε_xx[i] = z * (U[i-1] - 2*U[i] + U[i+1]) / abs(vozlisca[i+1,0] - vozlisca[i,0])*abs(vozlisca[i,0] - vozlisca[i-1,0])
    elif KE == 'TG1' or KE == 'TG2' or KE == 'TM': #Timošenko
        ε_xx = np.empty(N)
        ε_xx[0] = z * (-3/2*Φ[0] + 2*Φ[1] - 0.5*Φ[2]) / hx
        ε_xx[-1] = z * (3/2*Φ[-1] - 2*Φ[-2] + 0.5*Φ[-3]) / hx
        for i in range(1,N-1):
            ε_xx[i] = z * (-Φ[i-1] + Φ[i+1]) / (2*hx)
    return ε_xx
def epsilon_x_mke(z):
    # ε_xx = np.empty(N-1)
    # for i in range(N-1): ε_xx[i] = z/hx* (Φ[i+1] - Φ[i])
    ε_xx = z/hx* (Φ[1:] - Φ[:-1])
    return ε_xx
def epsilon_xz_mke(x):
    x = np.transpose([x])
    # ε_xz = 1/hx* (U[1:] - U[:-1] + hx*Φ[:-1] + x*(Φ[1:]-Φ[:-1])) # epsilon_xz = du/dx + phi NAPAKA!
    ε_xz = 1/hx* (U[1:] - U[:-1] - hx*Φ[:-1] - x*(Φ[1:]-Φ[:-1])) # epsilon_xz = du/dx - phi
    return ε_xz
# def gauss(U,m):
#     """Gaussova integracija"""
#     w = [[2], [1,1], [5/9,8/9,5/9], [0.347855,0.652145,0.347855,0.347855]] #uteži
#     x = [[0], [-1/np.sqrt(3),1/np.sqrt(3)], [-np.sqrt(3/5),0,np.sqrt(3/5)], [-0.861136,-0.339981,0.339981,0.861136]] #lege
#     integral = np.empty(N_e)
#     for el in range(N_e):
#         v1, v2 = vozlisca[elementi[el]]
#         I=0
#         for i in range(m):
#             u = U[el] + (U[el+1] - U[el]) / (v1[0] - v2[0]) * (x[m-1][i] - v1[0]) #linearna interpolacija
#             I += w[m-1][i]*u
#         integral[el] = I
#     return integral
# Ni potrebe po integraciji, ker imamo analitično rešitev - lin. potek / konst. po KE / debelini !

############## 1. deformacije
z = h/2 #največji odmik, kvadrat
# z = r #največji odmik, krog

# ε_xx_g = z * np.gradient(np.gradient(u[:,0], vozlisca[:,0], edge_order=2), vozlisca[:,0], edge_order=2) #slabše - napaka se sešteva
ε_xx_mkr = epsilon_x_mkr(z)
ε_xx = epsilon_x_mke(z)
X_KE = (vozlisca[1:,0] + vozlisca[:-1,0])/2 #središča KE

if KE == 'TG1' or KE == 'TG2' or KE == 'TM': #Timošenko
    γ_xz_mkr = np.empty(N)
    γ_xz_mkr[0] = -Φ[0] + (-3/2*U[0] + 2*U[1] - 0.5*U[2]) / hx
    γ_xz_mkr[-1] = -Φ[-1] + (3/2*U[-1] - 2*U[-2] + 0.5*U[-3]) / hx
    for i in range(1,N-1):
        γ_xz_mkr[i] = -Φ[i] + (-U[i-1] + U[i+1]) / (2*hx)

# γ_xz = epsilon_xz_mke([0,hx])
# γ_xz = np.ravel(γ_xz.T)
γ_xz = epsilon_xz_mke([hx/2])[0] #na sredini KE
X_KE2 = np.ravel([np.array([0,hx])+i*hx for i in range(N-1)]) #lin. potek po KE

ε_yy_mkr = ε_zz_mkr = -ν * ε_xx_mkr
ε_yy = ε_zz = -ν * ε_xx

γ_xz_analit = -F0/As/G*np.ones(50)

# plt.plot(vozlisca[:,0], np.gradient(u[:,0], vozlisca[:,0], edge_order=2), 'b.-')
plt.plot(vozlisca[:,0], ε_xx_mkr, 'k.--', label='$\\varepsilon_{xx}^{max}$ MKR')
plt.plot(X_KE, ε_xx, 'kx-', label='$\\varepsilon_{xx}^{max}$ MKE')
plt.plot(vozlisca[:,0], ε_yy_mkr, 'b.--', label='$\\varepsilon_{yy} = \\varepsilon_{zz}$ MKR')
plt.plot(X_KE, ε_yy, 'bx-', label='$\\varepsilon_{yy} = \\varepsilon_{zz}$ MKE')
if KE != 'EB': plt.plot(vozlisca[:,0], γ_xz_mkr, 'r.--', label='$\\gamma_{xz}$ MKR')
if KE != 'EB': plt.plot(X_KE, γ_xz, 'rx-', label='$\\gamma_{xz}$ MKE')
# if KE != 'EB': plt.plot(x, γ_xz_analit, 'C1--', label='$\\gamma_{xz}$ analit')
# plt.plot(vozlisca[:,0], ε_xx_g, 'r.-')
plt.xlabel('$x$ [mm]')
plt.ylabel('$\\varepsilon$ [/]')
plt.legend()
plt.minorticks_on()
plt.tight_layout()
plt.show()

############### 2. napetosti
σ_xx_mkr = E*ε_xx_mkr #elastično !!!
σ_xx = E*ε_xx #elastično !!!
σ_m_mkr = np.abs(σ_xx_mkr) #von Mises Euler-Bernoulli
σ_m = np.abs(σ_xx) #von Mises Euler-Bernoulli
if KE != 'EB': #Timošenko
    σ_xz_mkr = G*γ_xz_mkr
    σ_xz = G*γ_xz
    σ_m_mkr = np.sqrt(σ_xx_mkr**2 + 3*σ_xz_mkr**2) #von Mises
    σ_m = np.sqrt(σ_xx**2 + 3*σ_xz**2) #von Mises
    # σ_m = np.sqrt(np.repeat(σ_xx,2)**2 + 3*σ_xz**2) #von Mises
σ_xz_analit = -F0/As*np.ones(50) #iz analitične rešitve

plt.plot(vozlisca[:,0], σ_xx_mkr, 'k.--', label='$\sigma_{xx}^{max}$ MKR')
plt.plot(X_KE, σ_xx, 'kx-', label='$\sigma_{xx}^{max}$ MKE')
if KE != 'EB': 
    plt.plot(vozlisca[:,0], σ_xz_mkr, 'r.--', label='$\sigma_{xz}$ MKR')
    plt.plot(X_KE, σ_xz, 'rx-', label='$\sigma_{xz}$ MKE')
    # plt.plot(x, σ_xz_analit, 'C1-', label='$\sigma_{xz}$ analitično')
    plt.plot(vozlisca[:,0], σ_m_mkr, 'g.--', label='$\sigma_p$ von Mises MKR')
    plt.plot(X_KE, σ_m, 'gx-', label='$\sigma_p$ von Mises MKE')
    # plt.plot(vozlisca[:,0], σ_t, 'g.-', label='$\sigma_p$ Tresca')
if np.max(σ_m) > R_p: plt.axhline(R_p, c='gray', ls='-', label='$R_p$')
plt.xlabel('$x$ [mm]')
plt.ylabel('$\sigma$ [MPa]')
plt.legend()
plt.minorticks_on()
plt.tight_layout()
plt.show()

############# 3. notranje veličine
# dΦ_dx = np.gradient(Φ, vozlisca[:,0], edge_order=2) #OK
dΦ_dx = np.empty(N)
dΦ_dx[0] = (-3/2*Φ[0] + 2*Φ[1] - 0.5*Φ[2]) / hx
dΦ_dx[-1] = (3/2*Φ[-1] - 2*Φ[-2] + 0.5*Φ[-3]) / hx
for i in range(1,N-1):
    dΦ_dx[i] = (-Φ[i-1] + Φ[i+1]) / (2*hx)
# predpostavka sym prereza glede na nevtralno os
M_y = -E*I* dΦ_dx #notranji upogibni moment
if KE != 'EB': T_z = As*σ_xz_mkr #notranja strižna sila

# izpeljava iz KE [https://vulcanhammer.net/wp-content/uploads/2017/01/finite-elements-in-plasticity-theory-and-practice.pdf]
M_y_e = E*I/hx * (Φ[:-1] - Φ[1:]) #konst. po KE
if KE != 'EB': T_z_e = -G*As * ((U[1:] - U[:-1])/hx - (Φ[:-1] + Φ[1:])/2) #na sredini KE

# po definiciji MKE!
M_y_mke = -E*I/hx * (Φ[1:] - Φ[:-1]) #konst po KE
if KE != 'EB': T_z_mke = As*σ_xz #lin potek potek KE
# T_z_mke[:] = -F0*L #test analitična rešitev - niso prave notranje sile !!! TODO

# analitično za primer konzole levo, desno prosto + sila
x = np.linspace(0,L,100)
M_y_ex = F0*(L-x)
if KE != 'EB':T_z_ex = -F0*np.ones(len(x)) #TODO strig velika napaka!!!

fig, ax = plt.subplots()
ax.plot(vozlisca[:,0], M_y*1e-3, 'k.--', label='$M_y$ MKR')
# ax.plot(vozlisca[:-1,0], M_y_e*1e-3, 'kx', label='$M_y$ MKE')
ax.plot(X_KE, M_y_mke*1e-3, 'k+-', label='$M_y$ MKE')
ax.plot(x, M_y_ex*1e-3, 'k-.', label='analitično')
if KE != 'EB': #STRIG
    ax1=plt.twinx(ax)
    ax1.plot(vozlisca[:,0], T_z*1e-3, 'b.--',)
    # ax1.plot(vozlisca[:-1,0]+hx/2, T_z_e*1e-3, 'bx')
    ax1.plot(X_KE, T_z_mke*1e-3, 'b+-')
    ax1.plot(x, T_z_ex*1e-3, 'b-.', label='analitično')
    ax1.set_ylabel('$T_z$ [kN]', color='b')
    ax1.tick_params(axis='y', colors='b')
    alignYaxes([ax, ax1], [0,0])
ax.legend()
ax.set_xlabel('$x$ [mm]')
ax.set_ylabel('$M_y$ [Nm]')
plt.minorticks_on()
plt.tight_layout()
plt.show()

# globalni vektor vozliščnih notranjih obremenitev
f_int = np.zeros_like(F_glob)
for i in range(N-1):
    element = elementi[i]
    ind = np.array([Nps*element[0],Nps*element[0]+1,Nps*element[1],Nps*element[1]+1])
    # f_int[ind] += np.array([-T_z_mke[2*i], -M_y_mke[i], T_z_mke[2*i+1], M_y_mke[i]])
    f_int[ind] += np.array([-T_z_mke[i], -M_y_mke[i], T_z_mke[i], M_y_mke[i]]) #strig na sredini KE

r = F_glob - f_int #residual

# izločimo reakcije v podporah
if RP[0] == 1 or RP[0] == 2: #prvo vozlisce je fiksno
    F_Az = f_int[0]
    r[0] = 0
    print('F_Az = ', round(F_Az,2), 'N, analitično:', F0, 'N') #analitično le za KONZOLO!
if RP[0] == 2: #prvo vozlisce konzola
    M_Ay = f_int[1]
    r[1] = 0
    print('M_Ay = ', round(M_Ay*1e-3,2), 'Nm, analitično:', -F0*L*1e-3, 'Nm')
if RP[1] == 1 or RP[1] == 2: #zadnje vozlisce je fiksno
    F_Bz = f_int[-2]
    r[-2] = 0
    print('F_Bz = ', F_Bz)
if RP[1] == 2: #zadnje vozlisce konzola
    M_By = f_int[-1]
    r[-1] = 0
    print('M_By = ', M_By)

print('norm', np.linalg.norm(r), ', sum', np.sum(np.abs(r)), ', N = ',N)

fig, (ax1,ax2,ax3) = plt.subplots(3,1)
# s1=ax1.imshow(np.abs(F_ext).reshape(N,Nps).T, cmap='turbo')
s1=ax1.imshow(np.abs(F_glob).reshape(N,Nps).T, cmap='turbo')
s2=ax2.imshow(np.abs(f_int).reshape(N,Nps).T, cmap='turbo')
s3=ax3.imshow(np.abs(r).reshape(N,Nps).T, cmap='turbo')
plt.colorbar(s1,ax=ax1)
plt.colorbar(s2,ax=ax2)
plt.colorbar(s3,ax=ax3)
plt.grid(False)
plt.show()