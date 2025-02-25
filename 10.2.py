# -*- coding: utf-8 -*-
# 10.naloga: spektralno (direktno) reševanje poissonove enačbe za cilindrični KS v 2D

import numpy as np
import matplotlib.pyplot as plt
import scipy.fft as fft
from timeit import default_timer as timer

plt.style.use('aljaz')
plt.rcParams.update({'font.size': 14})
def pause():
    input("Press the <ENTER> key to continue...")

def TDMAsolver(a, b, c, d):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    and to http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
    '''
    nf = len(d) # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d)) # copy arrays
    for it in range(1, nf):
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1]
        dc[it] = dc[it] - mc*dc[it-1]

    xc = bc
    xc[-1] = dc[-1]/bc[-1]

    for il in range(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    return xc

R=1#enakostraničen valj 2R=H
H=2*R
N=50

#fiz prostor
z = np.linspace(0,H,N)
r = np.linspace(-R,R,N) #dz=dr
# r = np.linspace(-R,R,2*N) #dz=2*dr
r_j,z_k = np.meshgrid(r,z)
# j = np.arange(1,1+N,1,dtype=int)

#transformiran prostor
n=np.arange(1,1+N,1,dtype=int)
m=np.arange(1,1+N,1,dtype=int) #dz=dr
# m=np.arange(1,1+2*N,1,dtype=int)
m,n=np.meshgrid(m,n)

#RP
T1 = 10 #plašč, navodila T1>T2
T2 = -10 #zg,sp

t = timer()
# pas_sr = 2*np.cos(np.pi*n/N)-10 #če dz=2*dr
# pas_sp = 4 - 2/n
# pas_zg = 4 + 2/n

pas_sr = 2*np.cos(np.pi*n/N)-4 #če dz=dr
pas_sp = 1 - 1/2/n #n=j
pas_zg = 1 + 1/2/n
desna = np.zeros_like(pas_sr)

# RP r=-R
pas_sp[0] = 0
pas_sr[0] = 1
pas_zg[0] = 0
desna[0] = T1

# RP r=R
pas_sp[-1] = 0
pas_sr[-1] = 1
pas_zg[-1] = 0
desna[-1] = T1

desna = fft.dst(desna,axis=1) #delna transformacija (z-smer)
print(pas_sp.shape, pas_sr.shape, pas_zg.shape, desna.shape)

T_j_n = TDMAsolver(pas_sp.ravel()[:], pas_sr, pas_zg.ravel()[:], desna).reshape((N,N))
T_jk = fft.idst(T_j_n, axis=1) #delna transformacija (z-smer)
T_jk += T2 #Dirichlet RP v z smeri ????
print(T_jk.shape,T_j_n.shape)
t_komb = timer()-t
print(f'2.naloga: t = {t_komb:.3f} s')

# plot
pointsize=3
max_T = np.max(T_jk)

sc=plt.contourf(r_j,z_k,T_jk,cmap='coolwarm')
# # sc=plt.scatter(r_j,z_k,c=T_jk,s=pointsize, cmap='coolwarm')
cbar=plt.colorbar(sc,label='$T$')
# cbar_ticks = list(cbar.get_ticks()[:-1])  # Get existing ticks
# cbar_ticks.append(max_T)  # Add max(u) to the ticks
# cbar.set_ticks(cbar_ticks)

plt.axis('equal')
plt.grid(False)
plt.xlabel('$r$')
plt.ylabel('$z$')
plt.title(f'$N={N*N}$')
plt.tight_layout()
#plt.savefig('./temp.pdf')
plt.show()