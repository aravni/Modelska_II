# -*- coding: utf-8 -*-
# ravna plošča z rotacijo

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import scipy as sc
from matplotlib.cm import ScalarMappable

rc({"text.usetex": True,
    'family': 'serif', 'serif': ['Computer Modern']})
plt.rcParams.update({'font.size': 10})

U = 1 #potencial elektrode, RP dirichlet
l = 1 #širina panelov - RP območje
phi = 0*np.pi/4 # kot rotacije vseh panelov

def potencial(x,y,l=1):
    '''singularno v x=0, y=0!!!'''
    return 1/(2*np.pi) * (-l + y*(np.arctan((x+l/2)/y) - np.arctan((x-l/2)/y)) + (x+l/2)/2*np.log((x+l/2)**2 + y**2) - (x-l/2)/2*np.log((x-l/2)**2 + y**2))

def plosca(N,l=1,phi=0):
    dx = l/(N) #širina enega panela
    x_p = np.arange(-l/2+dx/2,l/2,dx) #središča panelov
    x_nodes = np.arange(-l/2,l/2+dx,dx)
    R = np.array([[np.cos(phi), np.sin(phi)],[-np.sin(phi), np.cos(phi)]])

    A = np.zeros((N,N))
    b = np.zeros((N))
    for i in range(N):
        for j in range(N):
            #xi,yi = x_p[i], 1e-10 #ni rotacije
            xi,yi = x_p[i] * np.array([np.cos(phi), np.sin(phi)])
            xj,yj = x_p[j] * np.array([np.cos(phi), np.sin(phi)])
            xi, yi = R @ np.array([xi-xj,yi-yj])
            if yi == 0.:
                yi = 1e-10
            A[i,j] = potencial(xi, yi, l=dx) # konst korak
        b[i] = U

    #sigma = np.linalg.solve(A,b) #ni ok, singularno!
    #sigma = np.linalg.lstsq(A,b)[0] #gostota naboja
    sigma = sc.linalg.lstsq(A,b)[0] #najbolje
    c = 1/2 * np.sum(sigma*dx)

    M=100
    meja=3
    x = np.linspace(-meja,meja,M)
    y = np.linspace(-meja,meja,M)
    
    u = np.zeros((M,M))
    for i in range(M):
        for j in range(M):
            for n in range(N):
                #xi, yi = x[i], y[j] #ni rotacije
                xi, yi = R @ np.array([x[i],y[j]])
                u[j,i] += potencial(xi-x_p[n], yi-0, l=dx) * sigma[n] #vrstica, stolpec

    # plotting
    plt.plot(x_p,sigma,'.-',label=f'{N}')
    plt.xlabel('elektroda')
    plt.ylabel('$\sigma$')

    # cbar=plt.contourf(x,y,u, vmin=-1, vmax=1)
    # #cbar=plt.contourf(x,y,u)
    # plt.axis('equal')
    # plt.title(f'$c=${c:.2f}, $N=${N}')
    # #plt.colorbar(cbar_vse, label='u')
    # plt.plot(x_nodes*np.cos(phi),x_nodes*np.sin(phi),'k.-')
    return c, sigma

NN = 3 #št izračunov po osi
N_list = np.linspace(1,NN*NN, NN*NN,dtype=int)
N_list = [1,2,3,5,10,20,40]

#fig, ax = plt.subplots(3,3, figsize=(12,10))
cbar_vse = ScalarMappable(cmap='viridis')
cbar_vse.set_array([])  # dummy mappable array
cbar_vse.set_clim(-1, 1)

c = []
for i in range(len(N_list)):
    print(f'{i+1}/{NN*NN}')
    #plt.subplot(NN,NN, int(1+i))
    ci, sigma = plosca(N=N_list[i], phi=phi)
    c.append(ci)
#fig.colorbar(cbar_vse, ax=ax, label='u')
plt.suptitle('različno število panelov')
plt.legend(title='$N$')
#plt.savefig('./sigma_plosca.pdf')
plt.show()

plt.plot(N_list, c, '.-')
plt.grid()
plt.xlabel('$N$')
plt.ylabel('c')
plt.title('odvisnost naboja na plošči od števila panelov')
#plt.savefig('./gost_naboja_paneli.pdf')
plt.show()

N=N_list[-1]#najveečje št panelov
plt.plot([dx*i for i in range(N)], sigma)
plt.show()

######################## doprinosi posameznih panelov
exit()
N = 4 #sodo št.
dx = l/(N) #širina enega panela
x_p = np.arange(-l/2+dx/2,l/2,dx) #središča panelov

A = np.zeros((N,N))
b = np.zeros((N))
for i in range(N):
    for j in range(N):
        #xi,yi = x_p[i], 1e-10 #ni rotacije
        xi,yi = x_p[i] * np.array([np.cos(phi), np.sin(phi)])
        if yi == 0.:
            yi = 1e-10
        A[i,j] = potencial(xi,yi, l=dx) # konst korak brez zasukov
    b[i] = U

sigma = sc.linalg.lstsq(A,b)[0]
c = 1/2 * np.sum(sigma*dx)

M=100
meja=5
x = np.linspace(-meja,meja,M)
y = np.linspace(-meja,meja,M)

R = np.array([[np.cos(phi), np.sin(phi)],[-np.sin(phi), np.cos(phi)]])
u = np.zeros((M,M,N))
for i in range(M):
    for j in range(M):
        for n in range(N):
            #xi, yi = x[i], y[j] #ni rotacije
            xi, yi = R @ np.array([x[i],y[j]])
            u[j,i,n] = potencial(xi-x_p[n],yi-0,l=dx) * sigma[n] #vrstica, stolpec!

# fig, ax = plt.subplots(N//2,N//2, figsize=(9,7))
# for n in range(N):
#     plt.subplot(N//2,N//2, int(1+n))
#     cbar = plt.contourf(x,y,u[:,:,n])
#     plt.axis('equal')
#     plt.title(f'panel {n+1}')
# fig.colorbar(cbar, ax=ax, label='u')
# plt.suptitle('potenciali posameznih panelov')
# plt.show()

fig, ax = plt.subplots()
for n in range(N):
    cbar = plt.contour(x,y,u[:,:,n], levels=10)
    plt.axis('equal')
fig.colorbar(cbar, ax=ax, label='u')
plt.suptitle(f'potenciali {N} panelov, $c=${c:.2f}')
#plt.savefig('./doprinos_panelov.pdf')
plt.show()