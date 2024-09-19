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

def potencial(x,y,l=1):
    '''singularno v y=0!!!'''
    return 1/(2*np.pi) * (-l + y*(np.arctan((x+l/2)/y) - np.arctan((x-l/2)/y)) + (x+l/2)/2*np.log((x+l/2)**2 + y**2) - (x-l/2)/2*np.log((x-l/2)**2 + y**2))
def R(phi):
    return np.array([[np.cos(phi), np.sin(phi)],[-np.sin(phi), np.cos(phi)]])

def interpolate_nodes(x_nodes, y_nodes, segments_per_line=1,closed=True):
    # Ensuring the array is closed
    if closed:
        x_nodes = np.append(x_nodes, x_nodes[0])
        y_nodes = np.append(y_nodes, y_nodes[0])
    
    # Prepare the interpolated arrays
    x_interp = []
    y_interp = []
    
    # Interpolating between nodes
    for i in range(len(x_nodes) - 1):
        x_line = np.linspace(x_nodes[i], x_nodes[i + 1], segments_per_line + 1)
        y_line = np.linspace(y_nodes[i], y_nodes[i + 1], segments_per_line + 1)
        
        # We skip the last point because it will be the first point of the next segment
        x_interp.extend(x_line[:-1])
        y_interp.extend(y_line[:-1])
        
    # Adding the last point of the last segment to close the loop
    #if closed:
    x_interp.append(x_nodes[-1])
    y_interp.append(y_nodes[-1])

    return np.array(x_interp), np.array(y_interp)

def plosca(x_nodes, y_nodes,closed=False):
    '''tvori območje z paneli na vozliščnih koordinatah'''
    N=len(x_nodes)
    if closed: start = 0
    else: start = 1

    # get middle points from nodes
    x_p = []
    y_p = []
    l_p = []
    phi_p = []
    for i in range(start,N):
        xi = (x_nodes[i-1]+x_nodes[i])/2
        yi = (y_nodes[i-1]+y_nodes[i])/2
        li = ((x_nodes[i-1]-x_nodes[i])**2 + (y_nodes[i-1]-y_nodes[i])**2)**0.5
        phi_i = np.arctan2(y_nodes[i]-y_nodes[i-1], x_nodes[i]-x_nodes[i-1])
        if yi == 0.: #singularno
            yi = 1e-10
        x_p.append(xi)
        y_p.append(yi)
        l_p.append(li)
        phi_p.append(phi_i)
    x_p = np.array(x_p)
    y_p = np.array(y_p)
    l_p = np.array(l_p)
    phi_p = np.array(phi_p)
    N=len(x_p) #št. panelov

    A = np.zeros((N,N))
    b = np.zeros((N))
    for i in range(N):
        for j in range(N):
            phi_i = np.arctan2(y_p[i]-y_p[j], x_p[i]-x_p[j]) #vrstni red i,j?!
            xi, yi = R(phi_i) @ np.array([x_p[i]-x_p[j],y_p[i]-y_p[j]])
            if yi == 0.:
                yi = 1e-10
            A[i,j] += potencial(xi, yi, l=l_p[i])
        b[i] = U
    sigma = np.linalg.solve(A,b) #ni ok, singularno!
    #sigma = np.linalg.lstsq(A,b)[0] #gostota naboja
    #sigma = sc.linalg.lstsq(A,b)[0] #najbolje
    c = 1/2 * np.sum(sigma*l_p)

    M=100 #to ja drago!
    meja=3
    x = np.linspace(-meja,meja,M)
    y = np.linspace(-meja,meja,M)

    u = np.zeros((M,M))
    for i in range(M):
        for j in range(M):
            for n in range(N):
                xi, yi = R(phi_p[n]) @ np.array([x[i]-x_p[n],y[j]-y_p[n]]) #globalni -> lokalni KS
                u[j,i] += potencial(xi, yi, l=l_p[n]) * sigma[n] #vrstica, stolpec

    # plotting
    # line_length = [np.sum(l_p[:i])+l_p[i]/2 for i in range(N)]
    # plt.plot(line_length, sigma,'.-', label=f'{N}')
    # plt.xlabel('dolžina elektrode')
    # plt.ylabel('$\sigma$')
    # #plt.title(f'$c=${c:.2f}, $N=${N}')
    # plt.xlim(0,np.sum(l_p))

    cbar=plt.contourf(x,y,u)
    plt.axis('equal')
    plt.title(f'$c=${c:.2f}, $N=${N}')
    plt.colorbar(cbar, label='u')
    #plt.plot(x_p,y_p,'ko')
    if closed:
        plt.plot(np.append(x_nodes, x_nodes[0]), np.append(y_nodes, y_nodes[0]), 'k.-')
    else:
        plt.plot(x_nodes,y_nodes, 'k.-')

############## 1. nodes kvadrat
# x0 = np.array([-1,1,1,-1])
# y0 = np.array([-1,-1,1,1])
# fig, ax = plt.subplots(1,3, figsize=(11,3))

# x,y = interpolate_nodes(x0, y0, 3)
# plt.subplot(1,3,1)
# plosca(x, y, False)

# x,y = interpolate_nodes(x0, y0, 10)
# plt.subplot(1,3,2)
# plosca(x, y, False)

# x,y = interpolate_nodes(x0, y0, 25)
# plt.subplot(1,3,3)
# plosca(x, y, False)

# plt.tight_layout()
# # plt.legend(title='$N$')
# #plt.savefig('./kvadrat.pdf')
# # plt.savefig('./kvadrat_sigma.pdf')
# plt.show()

# odprta stranica
# x = np.array([-1,1,1,-1])
# y = np.array([-1,-1,1,1])
# x,y = interpolate_nodes(x, y, 30, closed=False)
# plosca(x, y, closed=False)
# #plt.legend(title='$N$')
# plt.tight_layout()
# plt.savefig('./odprt_kvadrat.pdf')
# plt.show()

############## 2. nodes trikotnik
# x0 = np.array([1/2, -1/2, 0])
# y0 = np.array([-3**0.5/2/3, -3**0.5/2/3, 3**0.5/2*2/3])
# fig, ax = plt.subplots(1,3, figsize=(11,3))

# x,y = interpolate_nodes(x0, y0, 3)
# plt.subplot(1,3,1)
# plosca(x, y)

# x,y = interpolate_nodes(x0, y0, 10)
# plt.subplot(1,3,2)
# plosca(x, y)

# x,y = interpolate_nodes(x0, y0, 20)
# plt.subplot(1,3,3)
# plosca(x, y)

# plt.tight_layout()
# #plt.legend(title='$N$')
# plt.savefig('./trikot.pdf')
# plt.show()

############# 3. nodes krog
fig, ax = plt.subplots(1,3, figsize=(11,3))
N = 30
phi = np.linspace(0,2*np.pi,N,endpoint=True)
x = 1 * np.cos(phi)
y = 1 * np.sin(phi)

# očrtanje
# x2 = 0.5 * np.cos(phi)
# y2 = 0.5 * np.sin(phi)
# x = np.append(x,x2)
# y = np.append(y,y2)

plt.subplot(1,3,1)
plosca(x, y, closed=False)

N = 60
phi = np.linspace(0,2*np.pi,N,endpoint=True)
x = 1 * np.cos(phi)
y = 1 * np.sin(phi)
plt.subplot(1,3,2)
plosca(x, y, closed=False)

N = 100
phi = np.linspace(0,2*np.pi,N,endpoint=True)
x = 1 * np.cos(phi)
y = 1 * np.sin(phi)
plt.subplot(1,3,3)
plosca(x, y, closed=False)

plt.tight_layout()
#plt.legend(title='$N$')
#plt.savefig('./krog_dva.pdf')
plt.show()