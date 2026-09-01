#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Zaključna naloga: upogib nosilca, statika
# Euler-Bernoulli in Timošenkov nosilec, PLASTIČNOST
# viri:
# [1] E. A. de Souza Neto, D. Peric, D. R. J. Owen - Computational Methods for Plasticity: Theory and Applications - Wiley (2008)
# [2] D. R. J. Owen - Finite Elements in Plasticity: Theory and Practice - Pineridge Press Limited (1980) https://vulcanhammer.net/wp-content/uploads/2017/01/finite-elements-in-plasticity-theory-and-practice.pdf

import numpy as np 
import matplotlib.pyplot as plt
from scipy import linalg
from scipy.sparse.linalg import spsolve
from timeit import default_timer as timer
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

N = 10 # št. vozlišč, TODO zakaj samo sodo število!?
Nh = 1 # št. točk po višini prereza
KE = 'EB' #izbira KE (EB=nosilec Euler-Bernoulli, TG1,TG2,TM-oblike za Timošenko)

# parametri za Newton-Raphson-ovo metodo
tol = 1e-3
Nmax = 10

# snovne lastnosti
R_p = 300 # meja plastičnosti, MPa
ρ = 1 # kg/mm^3
E = 189126 # MPa
ν = 0.3 # Poisson
G = E/(2*(1+ν)) # strižni modul v elast. območju, MPa
Et = 0.003*E # konst. tangentni modul v plastičnem, MPa (np.nan za nelinearno)
H=H0 = Et/(1-Et/E) # linearno utrjevanje
#TODO modul utrjevanja H !!!

# obremenitve
RP = [2,0] # vpetje (0-prosto, 1-fiksno, 2-konzolno) [levo, desno]
RP_F = [0,0,1] # točkovna sila F0 (0-neobremenjeno, 1-obremenjeno) [levo, sredina, desno]
RP_M = [0,0,0] # upogibni moment M0 (0-neobremenjeno, 1-obremenjeno) [levo, sredina, desno]
F0 = 100000 # točkovna sila - navzdol, N
M0 = 500*1e3 # upogibni moment - smer ure, N mm
q0 = 0 # linijska obremenitev, N/mm

# geometrijske lastnosti
L = 1000 # dolžina nosilca, mm
# 1. pravokotnik: a-širina, h-višina
a, h = 100, 100 #mm
κ = 10*(1+ν)/(12+11*ν) #5/6 #koef. oblike prereza
A = a*h #mm^2
I = 1/12 * a*h**3 #mm^4
# 2. krog: r-polmer
# r = 100 #mm
# κ = 6*(1+ν)/(7+6*ν) #9/10 #koef. oblike prereza
# A = np.pi*r**2 #mm^2
# I = 1/64 * np.pi*r**4 #mm^4

hx = L/(N-1) # dolžina KE
As = κ*A # vpliv oblike prereza na velikost max strižnih nap., upoštevajoč konst. porazd. stižnih nap.
M_pl = R_p * a * h**2/4 # notranji moment pri meji plastičnosti

def remove_neg(a,b):
    # if a[0] <= 0:
    #     a = a[1:]
    #     b = b[1:]
    #     return remove_neg(a,b)
    # else: return a,b
    for i in range(len(a)): #remove ALL negative values
        if a[i] <= 0:
            a = np.delete(a, i)
            b = np.delete(b, i)
            return remove_neg(a,b)
    return a,b

# plastičnost - krivulja
εn = np.linspace(0,1,1000)
σn=[] #inž. napetost, nominalna
for i in range(len(εn)):
    σ_e = E*εn[i] #elastično
    if σ_e <= R_p: σn.append(σ_e)
    # else: σn.append(R_p) #idealno plastično -> Et=0
    else: σn.append(R_p + Et*(εn[i]-R_p/E)) #konstanten tangentni modul
σn = np.asarray(σn)

# # podatki iz LNMS
# # ROČNO pravi modul elastičnosti E in meja plastičnosti R_p!
# εn = np.load('./Modelska_II/zaklj_nal/epsilon.npy')
# σn = np.load('./Modelska_II/zaklj_nal/sigma.npy')

σt = σn*(1+εn) #Cauchy, prava napetost/true stress
ε_log = np.log(1+εn) #logaritemska deformacija
ε_pl = ε_log - σt/E #plastična deformacija
ε_pl, σt = remove_neg(ε_pl, σt)
    
# plt.plot(εn, σn, 'k-', label='$\sigma_{ing}$')
# plt.plot(ε_pl, σt, 'b-', label='$\sigma_{true}$')
# plt.legend()
# plt.xlabel('$\\varepsilon$ [/]')
# plt.ylabel('$\\sigma$ [MPa]')
# plt.show()

# reševanje sistema enačb za pomike
#############################################################################

# vozlišča in elementi
vozlisca = np.transpose([np.linspace(0,L,N), np.zeros(N)])
ind = np.linspace(0,N-1,N, dtype=int)
elementi = np.transpose([ind[:-1], ind[1:]]).astype(int)
Ne = len(elementi) #=N-1 št. KE

# sestav globalne masne in togostne matrike
Nps = 2
M_glob = np.zeros((N*Nps,N*Nps))
K_glob = np.zeros((N*Nps,N*Nps))
F_glob = np.zeros(N*Nps) # desna stran - vozliščne obremenitve

def distance(a,b): return ((a[0]-b[0])**2 + (a[1]-b[1])**2)**0.5
for element in elementi:
    Le = distance(vozlisca[element[0]], vozlisca[element[1]]) # dolžina elementa (Le = hx)
    # priševanje prispevkov posameznih elementov
    ind = np.array([[Nps*element[0],Nps*element[0]+1,Nps*element[1],Nps*element[1]+1]])
    if KE == 'EB': K_glob[ind.T,ind] += K_e_eb(E, I, Le)
    elif KE == 'TG1': K_glob[ind.T,ind] += K_e_tg1(E,I,Le)
    elif KE == 'TG2': K_glob[ind.T,ind] += K_e_tg2(E,I,Le)
    elif KE == 'TM': K_glob[ind.T,ind] += K_e_tm(E,I,Le)
    else: raise Exception('Napačen tip KE!')
    F_glob[ind[0,[0,2]]] += -q0*Le/2

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
t2=timer()
u = spsolve(K_glob, F_glob)
t3=timer()
U, Θ = u.reshape(vozlisca.shape).T
print(f'časi reševanja: sparse {(t3-t2):.3g} s')
print(f'maksimalni poves: {np.max(np.abs(U)):.5g} mm')

# # izris - elastično!
# # nedeformirano stanje
# [plt.plot(vozlisca[[i,j],0], vozlisca[[i,j],1],'-', c='gray') for i,j in elementi]
# plt.plot(vozlisca[:,0], vozlisca[:,1],'o', c='gray')
# # deformirano stanje
# [plt.plot(vozlisca[[i,j],0], U[[i,j]],'-', c='k') for i,j in elementi]
# plt.plot(vozlisca[:,0], U,'o', c='k')
# plt.axis('equal')
# plt.xlabel('$x$ [mm]')
# plt.ylabel('$y$ [mm]')
# plt.minorticks_on()
# plt.tight_layout()
# plt.show()

# plastičnost
###########################################################################################

def epsilon_x(z):
    '''z - višina v prerezu nosilca glede na nevtralno os'''
    ε_xx = z/hx* (Θ[1:] - Θ[:-1])
    return np.transpose([ε_xx])
def epsilon_xz(x):
    '''x - položaj po dolžini nosilca'''
    x = np.transpose([x])
    ε_xz = 1/hx* (U[1:] - U[:-1] - hx*Θ[:-1] - x*(Θ[1:]-Θ[:-1])) # epsilon_xz = du/dx - phi
    return ε_xz
def σ_y(ε_P):
    if Et is not np.nan:
        return R_p + H0*ε_P #linearno utrjevanje
    idx = np.argmin(np.abs(ε_pl - ε_P))
    σ_y = σt[idx]
    if σ_y < R_p:
        return R_p #elastično
    else:
        return σ_y
def hardening(ε_P):
    if Et is not np.nan:
        return H0 #linearno utrjevanje
    idx = np.argmin(np.abs(ε_pl - ε_P))
    σ_y = σt[idx]
    if σ_y < R_p:
        return E #elastično
    else:
        H = (σt[idx+1]- σt[idx-1]) / (ε_pl[idx+1] - ε_pl[idx-1])
        return H
def internal_f(σ, Θ):
    # notranje veličine
    M_y = -E*I/hx * (Θ[1:] - Θ[:-1]) #konst po KE
    if KE != 'EB': T_z = As*σ[1,:,0] #lin potek KE Timošenko
    else: T_z = np.zeros_like(M_y) #ni striga Euler-Bernoulli

    # globalni vektor notranjih vozliščnih sil
    f_int = np.zeros_like(F_glob)
    for i in range(N-1):
        element = elementi[i]
        ind = np.array([Nps*element[0],Nps*element[0]+1,Nps*element[1],Nps*element[1]+1])
        f_int[ind] += np.array([-T_z[i], -M_y[i], T_z[i], M_y[i]]) #strig na sredini KE
    r = F_glob - f_int #residual

    # izločimo reakcije v podporah TODO je to tukaj?
    if RP[0] == 1 or RP[0] == 2: #prvo vozlisce je fiksno
        F_Az = f_int[0]
        r[0] = 0
    if RP[0] == 2: #prvo vozlisce konzola
        M_Ay = f_int[1]
        r[1] = 0
    if RP[1] == 1 or RP[1] == 2: #zadnje vozlisce je fiksno
        F_Bz = f_int[-2]
        r[-2] = 0
    if RP[1] == 2: #zadnje vozlisce konzola
        M_By = f_int[-1]
        r[-1] = 0
    return f_int, r
def mises(sigma):
    if KE == 'EB': return np.abs(sigma[0])
    else: return np.sqrt(sigma[0]**2 + sigma[1]**2)
def local(ε_new, ε_old, εe_old, εe_new, εp_eq, σ):
    """lokalna zanka za von mises kriterij tečenja z linearnim izotropnim utrjevanjem"""
    # elastična napoved
    Δε = ε_new - ε_old #prirastek celotne deformacije je elastične (predpostavka)
    εp_eq_trial = εp_eq
    εe_trial = εe_old + Δε
    if KE == 'EB': σ_trial = E_mat*εe_trial #enoosno -> skalarno
    else: σ_trial = E_mat@εe_trial # Hooke
    σp_trial = mises(σ_trial) #mises primerjalna napetost
    Φ_trial = σp_trial - σ_y(εp_eq_trial)

    Φ = np.zeros_like(Φ_trial)
    Δγ = np.zeros_like(Φ_trial)
    conv=[]
    for i in range(Ne): #zanka po KE
        for j in range(Nh): #zanka po višini
            if Φ_trial[i,j] < 0: # elastično
                σ[:,i,j] = σ_trial[:,i,j]
                εe_new[:,i,j] = εe_trial[:,i,j]
                # εp = ε_old - εe_old #enako kot prejšnji korak
                εp_eq[i,j] = εp_eq_trial[i,j]
            else: # plastično
                l=0
                Φ[i,j] = Φ_trial[i,j]
                while abs(Φ[i,j]) > tol and l < Nmax:
                    # H = hardening(ε_P_i) #TODO nelin naklon!
                    Δγ[i,j] += Φ[i,j]/(3*G+H) #konstanten H za začetek
                    Φ[i,j] = σp_trial[i,j] - 3*G*Δγ[i,j] - σ_y(εp_eq_trial[i,j] + Δγ[i,j])
                    l += 1
                    conv.append(Φ[i,j])
                p = σ_trial[0,i,j] #hidrostatični del napetostnega tenzorja
                σd_trial = σ_trial[:,i,j] - p #TODO Euler Bernoulli? ni deviatorja!
                σd = (1 - Δγ[i,j]*3*G/σp_trial[i,j]) * σd_trial #deviatorični del napetostnega tenzorja
                εe_new[:,i,j] = 1/2/G * σd + 1/3*εe_trial[0,i,j] #elastični del deformacijskega tenzorja
                εp_eq[i,j] += Δγ[i,j]
                σ[:,i,j] = σd + p #napetostni tenzor
    return εe_new, εp_eq, σ, conv

# krivulja obremenjevanja
λ0 = 1e-6
N_incr = 10
λ = np.linspace(λ0,1,N_incr)

# točka izračuna (največjih) deformacij TODO več točk po višini prereza!!!
z = h/2 #največji odmik, kvadrat
# z = r #največji odmik, krog

# tenzor elastičnih konsatnt
if KE == 'EB': E_mat = np.array([E]) #enoosni Hooke
else: E_mat = np.array([[E,0],[0,G]]) #dodan strig, Hooke
# E_mat = np.array([[E,0],[0,G]]) #dodan strig, Hooke TODO v tenzor 4 reda ali Voigt??

# inicializacija
U = np.zeros(N) #začetni približek - lahko bolje!?
Θ = np.zeros(N)
f_int = np.zeros_like(F_glob) #uskladiti z neničelnim začetnim približkom - Hooke? TODO
# ε_new = np.zeros((3,3,Ne,Nh)) # TODO vse tenzorsko
if KE == 'EB': ε_new = np.zeros((1,Ne,Nh)) #dim {ε_xx}, EB ni striga
else: ε_new = np.zeros((2,Ne,Nh)) #dim {ε_xx, γ_xz}
εe_new = np.zeros_like(ε_new)
εp_eq = np.zeros((Ne,Nh)) #na začetku ni plastične deformacije
σ = np.zeros_like(ε_new) #napetostni tenzor

F_glob_norm = np.linalg.norm(F_glob)
conv_gl = []
conv_lok = []

for n in range(len(λ)): #koračno obremenjevanje
    print('loading state',n+1,'/',len(λ))
    f_ext = λ[n]*F_glob
    r = f_ext - f_int
    r_norm = np.linalg.norm(r) +F_glob_norm #da iteriramo!
    k=0
    while r_norm/F_glob_norm > tol and k < Nmax: #GLOBALNA ZANKA
        print('k=',k+1)
        # TODO posodobitev togosti K_glob!!!

        # Newton Raphson popravki
        δU, δΘ = spsolve(K_glob, r).reshape(vozlisca.shape).T
        U += δU
        Θ += δΘ
        
        # posodobitev deformacij
        ε_old = np.copy(ε_new)
        εe_old = np.copy(εe_new)
        ε_xx = epsilon_x(z) #TODO več točk po višini prereza
        if KE == 'EB': ε_new =  np.array([ε_xx])
        else:
            γ_xz = epsilon_xz([hx/2]) # sredina KE
            ε_new = np.array([ε_xx, γ_xz])

        εe_new, εp_eq, σ, conv = local(ε_new, ε_old, εe_old, εe_new, εp_eq, σ) # LOKALNA ZANKA

        f_int, r = internal_f(σ, Θ)
        r_norm=np.linalg.norm(r)
        conv_gl.append(abs(r_norm))
        conv_lok.extend(np.abs(conv))
        k += 1

conv_lok=np.array(conv_lok)
conv_gl=np.array(conv_gl)//np.linalg.norm(F_glob)
print(conv_gl.shape, conv_lok.shape)

# plt.plot(conv_gl, 'k.-', label='global')
plt.plot(conv_lok, 'b.-', label='lokal')
plt.yscale('log')
plt.xlabel('iteracija')
plt.ylabel('residual [/]')
plt.legend()
plt.tight_layout()
plt.show()

############## 1. deformacije
z = h/2 #največji odmik, kvadrat
# z = r #največji odmik, krog

ε_xx = epsilon_x(z)
# γ_xz = epsilon_xz([0,hx])
γ_xz = epsilon_xz([hx/2])[0] # sredina KE

# pretvorba v vozliščne vrednosti za vsak KE posebaj
ε_xx = np.repeat(ε_xx,2)
# γ_xz = np.ravel(γ_xz.T)
γ_xz = np.repeat(γ_xz,2)
ε_yy = ε_zz = -ν * ε_xx
ε_E = np.transpose([ε_xx, γ_xz])

X_KE = (vozlisca[1:,0] + vozlisca[:-1,0])/2 #središča KE
X_KE2 = np.ravel([np.array([0,hx])+i*hx for i in range(N-1)]) #lin. potek po KE

# plt.plot(vozlisca[:,0], np.gradient(u[:,0], vozlisca[:,0], edge_order=2), 'b.-')
plt.plot(X_KE2, ε_xx, 'k.-', label='$\\varepsilon_{xx}^{max}$')
plt.plot(X_KE2, ε_yy, 'b.-', label='$\\varepsilon_{yy} = \\varepsilon_{zz}$')
if KE != 'EB': plt.plot(X_KE2, γ_xz, 'r.-', label='$\\gamma_{xz}$')
# plt.plot(vozlisca[:,0], ε_xx_g, 'r.-')
plt.xlabel('$x$ [mm]')
plt.ylabel('$\\varepsilon$ [/]')
plt.legend()
plt.minorticks_on()
plt.tight_layout()
plt.show()

############### 2. napetosti
σ_xx = E*ε_xx #elastično !!!
σ_m0 = np.abs(σ_xx) #von Mises Euler-Bernoulli
σ_xz = np.zeros_like(σ_xx)
if KE != 'EB': #Timošenko
    σ_xz = G*γ_xz
    σ_m0 = np.sqrt(σ_xx**2 + 3*σ_xz**2) #von Mises
σ_E = np.transpose([σ_xx, σ_xz])

# popravek pomikov vsak inkrement, glede na residual sile na desni strani r=f_int-f_ext
conv=[]

for n in range(len(λ)):
    print('loading state',n+1,'/',len(λ))
    # popravki neznank
    δU, δΘ = spsolve(K_glob, r).reshape(vozlisca.shape).T
    U += δU
    Θ += δΘ

    # deformacije
    ε_xx = epsilon_x(z)
    γ_xz = epsilon_xz([hx/2])[0] # sredina KE

    # pretvorba v vozliščne vrednosti za vsak KE posebaj
    ε_xx = np.repeat(ε_xx,2)
    γ_xz = np.repeat(γ_xz,2)
    ε_yy = ε_zz = -ν * ε_xx
    ε_E = np.transpose([ε_xx, γ_xz])

    # napetosti
    σ_xx = E*ε_xx #elastično !!!
    σ_xz = G*γ_xz
    σ_m0 = np.sqrt(σ_xx,2**2 + 3*σ_xz**2) #von Mises
    σ_E = np.transpose([σ_xx, σ_xz])

    # notranje veličine
    M_y = -E*I/hx * (Θ[1:] - Θ[:-1]) #konst po KE
    T_z = As*σ_xz #lin potek potek KE

    # lokalna zanka
    σ = np.zeros((N,3,3))
    ε = np.zeros((N,3,3))
    ε_P = np.zeros((N))
    for i in range(N):
        # deformacijski tenzor
        ε_E_trial = np.array([[ε_E[i,0],    0,          ε_E[i,1]],
                            [0,           ε_yy[i],    0],
                            [ε_E[i,1],    0,          ε_zz[i]]])
        ε_P_trial = 0.0
        ε_v_trial = 1/3*(ε_xx[i] + ε_yy[i] + ε_zz[i]) #volumetrični/krogelni/hidrostatični del deformacijskega tenzorja
        ε_d_trial = ε_E_trial - ε_v_trial*np.eye(3) #deviatorični del deformacijskega tenzorja
        
        # elastično
        p_trial = E*ε_v_trial #hidrostatični
        s_trial = 2*G*ε_d_trial #deviatorični
        σ_trial = s_trial + p_trial*np.eye(3) #napetostni tenzor
        q_trial = np.sqrt(3/2 * np.tensordot(s_trial, s_trial))

        # test plastičnosti
        Φ_trial = q_trial - σ_y(ε_P_trial) #von Mises

        if Φ_trial < 0: #elastično
            # print(i+1, 'elastično')
            σ[i] = σ_trial
            ε[i] = ε_E_trial
        else: #plastično
            # print(i+1, 'plastično')
            
            Δγ = 0 #zač. pribl.

            ε_P_i = ε_P_trial
            s = s_trial
            q = q_trial
            Φ = Φ_trial

            k=0
            while np.sum(np.abs(Φ)) >= tol and k < Nmax:
                H = hardening(ε_P_i)
                
                Δγ += Φ/(3*G+H)
                ε_P_i += Δγ
                s = (1-Δγ*3*G/q)*s
                q = np.sqrt(3/2 * np.tensordot(s, s))
                Φ = q - σ_y(ε_P_i) #von Mises

                k+=1

            ε_P[i] = ε_P_i
            σ[i] = s + p_trial*np.eye(3) #napetostni tenzor
            ε[i] = 1/(2*G)*s + 1/3*ε_v_trial*np.eye(3) #deformacijski tenzor
        
    # notranje veličine
    M_y = -E*I/hx * (Θ[1:] - Θ[:-1]) #konst po KE
    T_z = As*σ[:,0,2] #lin potek potek KE

    # globalni vektor notranjih vozliščnih sil
    f_int = np.zeros_like(F_glob)
    for i in range(N-1):
        element = elementi[i]
        ind = np.array([Nps*element[0],Nps*element[0]+1,Nps*element[1],Nps*element[1]+1])
        f_int[ind] += np.array([-T_z[i], -M_y[i], T_z[i], M_y[i]]) #strig na sredini KE
    r = F_glob - f_int #residual

    # izločimo reakcije v podporah
    if RP[0] == 1 or RP[0] == 2: #prvo vozlisce je fiksno
        F_Az = f_int[0]
        r[0] = 0
    if RP[0] == 2: #prvo vozlisce konzola
        M_Ay = f_int[1]
        r[1] = 0
    if RP[1] == 1 or RP[1] == 2: #zadnje vozlisce je fiksno
        F_Bz = f_int[-2]
        r[-2] = 0
    if RP[1] == 2: #zadnje vozlisce konzola
        M_By = f_int[-1]
        r[-1] = 0

    # nova obremenitev
    f_ext = λ[n]*F_glob
    r = f_ext-f_int #residual
    conv.append(np.sum(np.abs(r)))

    σ_m = np.sqrt(σ[:,0,0]**2 + 3*σ[:,0,2]**2) #von Mises
    print('residual:', np.sum(np.abs(r)))
plt.semilogy(conv, 'k.-')
plt.show()
input('end')

σ = np.zeros((N,3,3))
ε = np.zeros((N,3,3))
ε_P = np.zeros((N))
conv=[]
for i in range(N):
    # deformacijski tenzor
    ε_E_trial = np.array([[ε_E[i,0],    0,          ε_E[i,1]],
                          [0,           ε_yy[i],    0],
                          [ε_E[i,1],    0,          ε_zz[i]]])
    ε_P_trial = 0.0
    ε_v_trial = 1/3*(ε_xx[i] + ε_yy[i] + ε_zz[i]) #volumetrični/krogelni/hidrostatični del deformacijskega tenzorja
    ε_d_trial = ε_E_trial - ε_v_trial*np.eye(3) #deviatorični del deformacijskega tenzorja
    
    # elastično
    p_trial = E*ε_v_trial #hidrostatični
    s_trial = 2*G*ε_d_trial #deviatorični
    σ_trial = s_trial + p_trial*np.eye(3) #napetostni tenzor
    q_trial = np.sqrt(3/2 * np.tensordot(s_trial, s_trial))

    # test plastičnosti
    Φ_trial = q_trial - σ_y(ε_P_trial) #von Mises

    if Φ_trial < 0: #elastično
        print(i+1, 'elastično')
        σ[i] = σ_trial
        ε[i] = ε_E_trial
    else: #plastično
        print(i+1, 'plastično')
        
        Δγ = 0 #zač. pribl.

        ε_P_i = ε_P_trial
        s = s_trial
        q = q_trial
        Φ = Φ_trial

        k=0
        convi=[np.abs(Φ)]
        while np.sum(np.abs(Φ)) >= tol and k < Nmax:
            H = hardening(ε_P_i)
            
            Δγ += Φ/(3*G+H)
            ε_P_i += Δγ
            s = (1-Δγ*3*G/q)*s
            q = np.sqrt(3/2 * np.tensordot(s, s))
            Φ = q - σ_y(ε_P_i) #von Mises

            k+=1
            convi.append(np.abs(Φ))
        conv.append(convi)

        ε_P[i] = ε_P_i
        # s = (1-Δγ*3*G/q_trial)*s_trial
        σ[i] = s + p_trial*np.eye(3) #napetostni tenzor
        ε[i] = 1/(2*G)*s + 1/3*ε_v_trial*np.eye(3) #deformacijski tenzor
σ_m = np.sqrt(σ[:,0,0]**2 + 3*σ[:,0,2]**2) #von Mises

f_int = np.zeros(N*Nps)
for i in range(Ne):
    element = elementi[i]
    ind = np.array([Nps*element[0],Nps*element[0]+1,Nps*element[1],Nps*element[1]+1])
    if i > 0 and i < Ne-1: #2. reda natančnosti
        f1 = A*G/2 * (Θ[i]+Θ[i+1] + (U[i+2]-U[i] + U[i+1]-U[i-1])/2/hx)
        f2 = a*h**2*E/8/hx * (Θ[i+1]-Θ[i-1] + Θ[i+2]-Θ[i]) #TODO mora biti tukaj E*I ?!
    elif i==0:
        f1 = A*G/2 * (Θ[i]+Θ[i+1] + (U[i+2]-U[i])/2/hx + (-0.5*U[i+2]+2*U[i+1]-3/2*U[i])/hx)
        f2 = a*h**2*E/4/hx * (-0.5*Θ[i+2]+2*Θ[i+1]-3/2*Θ[i] + (Θ[i+2]-Θ[i])/2)
    elif i==Ne-1:
        f1 = A*G/2 * (Θ[i]+Θ[i+1] + (3/2*U[i+1]-2*U[i]+0.5*U[i-1])/hx + (U[i+1]-U[i-1])/2/hx)
        f2 = a*h**2*E/4/hx * ((Θ[i+1]-Θ[i-1])/2 + 3/2*Θ[i+1]-2*Θ[i]+0.5*Θ[i-1])
    f_int[ind] += np.array([-f1,f2,f1,-f2])

print(f_int.shape, F_glob.shape, np.sum(np.abs(F_glob-f_int)))
plt.imshow(np.abs(F_glob-f_int).reshape(N,Nps).T, cmap='turbo')
plt.colorbar()
plt.title('napake v notranjih silah')
plt.xlabel('vozlišče')
plt.ylabel('komponenta')
plt.minorticks_on()
plt.tight_layout()
plt.show()

for i in range(len(conv)):
    plt.plot(conv[i], '.-')
plt.yscale('log')
plt.xlabel('iteracija')
plt.ylabel('residual')
plt.minorticks_on()
plt.tight_layout()
plt.show()

plt.plot(X_KE2, σ_xx, 'k.-', label='$\sigma_{xx}^{max}$')
plt.plot(vozlisca[:,0], σ[:,0,0], 'k.--', label='$\sigma_{xx}^n$')
if KE != 'EB': 
    plt.plot(X_KE2, σ_xz, 'r.-', label='$\sigma_{xz}$')
    plt.plot(vozlisca[:,0], σ[:,0,2], 'r.--', label='$\sigma_{xz}^n$')
    plt.plot(X_KE2, σ_m0, 'g.-', label='$\sigma_p$ von Mises')
    plt.plot(vozlisca[:,0], σ_m, 'g.--', label='$\sigma_p^n$ von Mises')
    # plt.plot(vozlisca[:,0], σ_t, 'g.-', label='$\sigma_p$ Tresca')
plt.axhline(R_p, c='gray', ls='-', label='$R_p$')
plt.xlabel('$x$ [mm]')
plt.ylabel('$\sigma$ [MPa]')
plt.legend()
plt.minorticks_on()
plt.tight_layout()
plt.show()

    # test plastičnosti
    # ξ_trial = 2*σ_trial[:,1]**2 + 2/3*σ_trial[:,0]**2
    # Φ_trial = ξ_trial/2 - 1/3*σ_y(ε_P_trial)**2

    # conv=[]
    # for i in range(N):
    #     if Φ_trial[i] < 0: #elastično
    #         print(i,'elastično')
    #         σ = σ_trial
    #         ε = ε_E_trial
    #     else: #plastično
    #         print(i,'plastično')
            
    #         tol=1e-5
    #         Nmax=1000
    #         Δγ = 0 #zač. pribl.

    #         ξ = 2/3*σ_trial[i,0]**2 + 2*σ_trial[i,1]**2
    #         Φ = ξ/2 - 1/3*σ_y(ε_P_trial)**2

    #         σ = σ_trial
    #         ε = ε_E_trial
    #         ε_P = ε_P_trial

    #         conv_i=[]
    #         k=0
    #         while np.abs(Φ) > tol and k < Nmax:
    #             ξ_ = -8*G*σ[i,1]**2 / (1+2*G*Δγ)**3 - 8/9*E*σ[i,0]**2 / (1+2/3*E*Δγ)**3
    #             H = Et #modul utrjevanja TODO nelin
    #             H_ = 2*σ_y(ε_P) * (2/3)**0.5*H * (ξ**0.5 + Δγ*ξ_/(2*ξ**0.5))
    #             Φ_= ξ_/2 - 1/3*H_

    #             Δγ += - Φ/Φ_ #popravek
    #             ε_P += Δγ*(2/3*ξ)**0.5 #popravek

    #             ξ = 2/3*σ_trial[i,0]**2 + 2*σ_trial[i,1]**2
    #             Φ = ξ/2 - 1/3*σ_y(ε_P)**2

    #             k+=1
    #             conv_i.append(np.abs(Φ))
    #             # print(Φ)
    #         conv.append(conv_i)
    #         # input(k)
    # for i in range(len(conv)):
    #     plt.plot(conv[i], '.-')
    # plt.yscale('log')
    # plt.xlabel('iteracija')
    # plt.ylabel('residual')
    # plt.minorticks_on()
    # plt.tight_layout()
    # plt.show()

    # if Φ_trial.all() < 0: #elastično
    #     print('elastično')
    #     σ = σ_trial
    #     ε = ε_E_trial
    # else: #plastično
    #     print('plastično')
        
    #     tol=1e-2
    #     Nmax=1000
    #     Δγ = 0 #zač. pribl.

    #     ξ = 2/3*σ_trial[:,0]**2 + 2*σ_trial[:,1]**2 #TODO FEM drugače!!!
    #     Φ = ξ/2 - 1/3*σ_y(ε_P_trial)**2

    #     σ = σ_trial
    #     ε = ε_E_trial
    #     ε_P = ε_P_trial

    #     k=0
    #     conv=[]
    #     while np.sum(np.abs(Φ)) >= tol and k < Nmax:
    #         ξ_ = -8*G*σ[:,1]**2 / (1+2*G*Δγ)**3 - 8/9*E*σ[:,0]**2 / (1+2/3*E*Δγ)**3
    #         H = Et #modul utrjevanja TODO nelin
    #         H_ = 2*σ_y(ε_P) * (2/3)**0.5*H * (ξ**0.5 + Δγ*ξ_/(2*ξ**0.5))
    #         Φ_= ξ_/2 - 1/3*H_

    #         ξ = 2/3*σ_trial[:,0]**2 + 2*σ_trial[:,1]**2
    #         Φ = ξ/2 - 1/3*σ_y(ε_P)**2

    #         Δγ += - Φ/Φ_ #popravek
    #         ε_P += Δγ*(2/3*ξ)**0.5 #popravek

    #         k+=1
    #         conv.append(np.sum(np.abs(Φ)))

############# 3. notranje veličine
M_y = -E*I/hx * (Θ[1:] - Θ[:-1]) #konst po KE
T_z = As*σ_xz #lin potek potek KE

# analitično za primer konzole levo, desno prosto + sila
x = np.linspace(0,L,100)
M_y_ex = F0*(L-x)
T_z_ex = -F0*L*np.ones(len(x))

print(M_y.shape, T_z.shape, F_glob.shape)
f_int = np.zeros(N*Nps)
for i in range(N-1):
    element = elementi[i]
    ind = np.array([Nps*element[0],Nps*element[0]+1,Nps*element[1],Nps*element[1]+1])
    f_int[ind] += np.array([-T_z[2*i], -M_y[i], T_z[2*i+1], M_y[i]])*1e-3
print(f_int.shape, F_glob.shape, np.sum(np.abs(F_glob-f_int)))
fig, (ax1,ax2,ax3) = plt.subplots(3,1)
s1=ax1.imshow(np.abs(F_glob).reshape(N,Nps).T, cmap='turbo')
s2=ax2.imshow(np.abs(f_int).reshape(N,Nps).T, cmap='turbo')
s3=ax3.imshow(np.abs(f_int-F_glob).reshape(N,Nps).T, cmap='turbo')
plt.colorbar(s1,ax=ax1)
plt.colorbar(s2,ax=ax2)
plt.colorbar(s3,ax=ax3)
plt.grid(False)
plt.show()

fig, ax = plt.subplots()
ax.plot(X_KE, M_y*1e-3, 'kx-', label='MKE')
ax.plot(x, M_y_ex*1e-3, 'k--', label='analitično')
# ax.axhline(M_pl*1e-3, c='gray', ls='-', label='$M_{pl}$')
ax1=plt.twinx(ax)
ax1.plot(X_KE2, T_z*1e-3, 'bx-')
ax1.plot(x, T_z_ex*1e-3, 'b--', label='analitično')
ax.legend(loc='upper left')
ax.set_xlabel('$x$ [mm]')
ax.set_ylabel('$M_y$ [Nm]')
ax1.set_ylabel('$T_z$ [kN]', color='b')
ax1.tick_params(axis='y', colors='b')
alignYaxes([ax, ax1], [0,0])
plt.minorticks_on()
plt.tight_layout()
plt.show()
