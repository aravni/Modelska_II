# -*- coding: utf-8 -*-
# spektralno reševanje poissonove enačbe za cilindrični KS v 2D

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import scipy.fft as fft
from timeit import default_timer as timer

rc({"text.usetex": True,
    'family': 'serif', 'serif': ['Computer Modern']})
plt.rcParams.update({'font.size': 10})
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

R=1#enakostraničen valj 2r=h
H=2*R
N=10

#fiz prostor
z = np.linspace(0,H,N) #dz=dr
r = np.linspace(-R,R,N)
r_j,z_k = np.meshgrid(r,z)
j = np.arange(1,1+N,1,dtype=int)

#transformiran prostor
m=n=np.arange(1,1+N,1,dtype=int)
m,n=np.meshgrid(m,n)

t = timer()
# pas_sr = 2*np.cos(np.pi*n/N)-10 #če dz=2*dr
# pas_sp = 4 - 2/j[1:]
# pas_zg = 4 + 2/j[:-1]

pas_sr = 2*np.cos(np.pi*n/N)-10 #če dz=dr
pas_sp = 4 - 2/n
pas_zg = 4 + 2/n
desna = np.zeros_like(pas_sr)

#RP
T1 = 0 #zg,sp
T2 = 1 #plašč, navodila T2>T1

i=0 #sp
pas_sr[i,:] = 1
pas_sp[i,:] = 0
pas_zg[i,:] = 0
desna[i,:] = T1

i=-1 #zg
pas_sr[i,:] = 1
pas_sp[i,:] = 0
pas_zg[i,:] = 0
desna[i,:] = T1

i=0 #
pas_sr[:,i] = 1
pas_sp[:,i] = 0
pas_zg[:,i] = 0
desna[:,i] = T2

i=-1 #
pas_sr[:,i] = 1
pas_sp[:,i] = 0
pas_zg[:,i] = 0
desna[:,i] = T2

desna = fft.dst(desna,axis=0) #delna transformacija

T_j_n = TDMAsolver(pas_sp.flatten()[:], pas_sr, pas_zg.flatten()[:], desna).reshape((N,N))
T_jk = fft.idst(T_j_n,axis=0) #delna transformacija
t_komb = timer()-t
print('kombinirana: t = ',t_komb)

# plot
pointsize=3
max_T = np.max(T_jk)

sc=plt.contourf(r_j,z_k,T_jk,cmap='coolwarm')
#sc=plt.scatter(x,y,c=u,s=pointsize)
cbar=plt.colorbar(sc,label='T')
cbar_ticks = list(cbar.get_ticks()[:-1])  # Get existing ticks
cbar_ticks.append(max_T)  # Add max(u) to the ticks
cbar.set_ticks(cbar_ticks)

plt.axis('equal')
plt.xlabel('$r$')
plt.ylabel('$z$')
plt.title(f'$N={T_jk.size}$')
plt.tight_layout()
#plt.savefig('./temp.pdf')
plt.show()