#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 8.naloga: MKE - lastne rešitve

import numpy as np 
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from scipy.linalg import eig, eigh
from scipy.sparse.linalg import eigs, eigsh
from timeit import default_timer as timer
import matplotlib
from scipy.spatial import KDTree
from scipy.special import jn, jn_zeros
from MESH import *
from nal8_galerkin import calc

matplotlib.use('TkAgg') # Qt čudne težave z Wayland
plt.style.use('aljaz')

def analit(m,n, r=None, phi=None):
    '''analitična rešitev za lastne oblike polkrožne opne, m,n >= 1'''
    if r is None and phi is None:
        r = np.linspace(0,1,100)
        phi = np.linspace(0,np.pi,180)
        r,phi = np.meshgrid(r,phi)
    x=r*np.cos(phi)
    y=r*np.sin(phi)
    j_mn=jn_zeros(m,n)[-1]
    u = jn(m,j_mn*r) * np.sin(m*phi)
    return x,y,u, j_mn
# x,y,u,k=analit(1,1)
# print('lastna vrednost k = ',k)
# plt.tricontourf(x.ravel(),y.ravel(),u.ravel(), cmap='turbo')
# plt.show()
def analit_sq(m,n, X,Y):
    '''analitična rešitev za lastne oblike kvadratne opne, m,n >= 1'''
    u = np.sin(m*np.pi*X) * np.sin(n*np.pi*Y) #X,Y meshgrid!
    k = np.sqrt((m*np.pi)**2 + (n*np.pi)**2)
    return X,Y,u,k

def solver(N, opna, mreza, Nc, Nr, plot=False):
    t1=timer()
    generate_mesh = True
    if opna == 'kvadrat': # kvadratna enotska opna (L=H=1)
        h=1/(N-1)
        x=y=np.linspace(0,1,N)
        X,Y=np.meshgrid(x,y)
        pointsall = np.array([X.flatten(), Y.flatten()]).T

        i=np.linspace(0,N-1,N, dtype=int)
        bc = np.unique([i, i+(N-1)*N, i*N, i*N+N-1]).astype(int) #indeksi robnih vozlišč

        figsize=(Nc*2+0.5, Nr*2)

    elif opna == 'f': # f opna 6.nal
        h=1/(N-1)
        x=y=np.linspace(0,1,N)
        X,Y=np.meshgrid(x,y)

        t = np.empty_like(X, dtype=bool) # TODO remove mesh from outside domain!! -> import mesh!!!
        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                if 1/4 < X[i,j] < 1/2:
                    t[i,j] = True
                elif Y[i,j] > 3*1/4 and X[i,j] > 1/4:
                    t[i,j] = True
                elif 1/4 < Y[i,j] < 1/2 and X[i,j] < 3/4:
                    t[i,j] = True
                else: t[i,j] = False
        X_out, Y_out = X[~t].ravel(), Y[~t].ravel()
        X=X[t]
        Y=Y[t]
        pointsall = np.array([X.flatten(), Y.flatten()]).T
        # bc=[]
        # for i in range(N): #PROBLEM: le zunanje meje!
        #     hi = h*i
        #     x_min = y_min = np.inf
        #     x_max = y_max = -np.inf
        #     for j in range(pointsall.shape[0]):
        #         if pointsall[j,1] == hi:
        #             if pointsall[j,0] > x_max:
        #                 x_max = pointsall[j,0]
        #                 i_x_max = j
        #             if pointsall[j,0] < x_min:
        #                 x_min = pointsall[j,0]
        #                 i_x_min = j
        #         if pointsall[j,0] == hi:
        #             if pointsall[j,1] > y_max:
        #                 y_max = pointsall[j,1]
        #                 i_y_max = j
        #             if pointsall[j,1] < y_min:
        #                 y_min = pointsall[j,1]
        #                 i_y_min = j
        #     bc.append(i_x_max)
        #     bc.append(i_x_min)
        #     bc.append(i_y_max)
        #     bc.append(i_y_min)

        tol = 1e-8
        tree = KDTree(pointsall) #robust way for finding if any neigbour is missing
        bc = []
        for i, (x, y) in enumerate(pointsall):
            neighbors = tree.query_ball_point([(x + h, y), (x - h, y), (x, y + h), (x, y - h)], r=tol)
            if any(len(n) == 0 for n in neighbors):
                bc.append(i)
        bc = np.array(bc).astype(int)

        figsize=(Nc*2+0.5, Nr*2)

    elif opna == 'polkrog': # polkrožna enotska opna (R=1)
        dr=1/(N-1)
        d_phi=2*np.pi/(N-1)
        if mreza == 'random polar': # naključna vozlišča
            R=np.random.rand(N,N)
            PHI=np.random.rand(N,N)*np.pi
            i=R.size
            R=np.append(R,0)
            PHI=np.append(PHI,0)
            R=np.append(R,np.full(N,1))
            PHI=np.append(PHI,np.linspace(0,np.pi,N))
            R=np.append(R,np.linspace(dr,1,N))
            PHI=np.append(PHI,np.full(N,0))
            R=np.append(R,np.linspace(dr,1,N))
            PHI=np.append(PHI,np.full(N,np.pi))
            X=R*np.cos(PHI)
            Y=R*np.sin(PHI)
            pointsall = np.array([X.flatten(), Y.flatten()]).T
            bc = np.linspace(i, R.size-1, R.size, dtype=int) #robna vozlišča
        elif mreza == 'polar':
            r=np.linspace(dr,1,N)
            phi=np.linspace(0,np.pi,N)
            R,PHI=np.meshgrid(r,phi)
            R=np.append(R,0) #le 1x izhodišče!
            PHI=np.append(PHI,0)
            X=R*np.cos(PHI)
            Y=R*np.sin(PHI)
            pointsall = np.array([X.flatten(), Y.flatten()]).T
            i=np.linspace(0,N-1,N, dtype=int)
            bc = np.unique([i, i+(N-1)*N, np.full_like(i,X.ravel().size-1), i*N+N-1]).astype(int)
        elif mreza == 'cartesian':
            x = np.linspace(-1+dr, 1-dr, N)
            y = np.linspace(dr, 1, N)
            X, Y = np.meshgrid(x, y)
            pointsall = np.vstack([X.ravel(), Y.ravel()]).T
            pointsall = pointsall[X.ravel()**2 + Y.ravel()**2 < 1]
            phi_boundary = np.linspace(0, np.pi, N)
            boundary_arc = np.column_stack((np.cos(phi_boundary), np.sin(phi_boundary)))
            diameter = np.column_stack((np.linspace(-1+dr, 1-dr, N), np.zeros(N)))
            i=pointsall.shape[0]
            pointsall = np.vstack([pointsall, boundary_arc, diameter]) #cartesian!
            X = pointsall[:, 0]
            Y = pointsall[:, 1]
            R = np.sqrt(X**2 + Y**2)
            PHI = np.arctan2(Y, X)
            bc = np.linspace(i, pointsall.shape[0]-1, pointsall.shape[0], dtype=int)
        elif mreza == 'random cartesian':
            x = np.random.rand(N, N) * 2 - 1
            y = np.random.rand(N, N)
            pointsall = np.vstack([x.ravel(), y.ravel()]).T
            pointsall = pointsall[np.sqrt(x.ravel()**2 + y.ravel()**2) < 1]
            phi_boundary = np.linspace(0, np.pi, N)
            boundary_arc = np.column_stack((np.cos(phi_boundary), np.sin(phi_boundary)))
            diameter = np.column_stack((np.linspace(-1+dr, 1-dr, N), np.zeros(N)))
            i=pointsall.shape[0]
            pointsall = np.vstack([pointsall, boundary_arc, diameter])
            X = pointsall[:, 0]
            Y = pointsall[:, 1]
            R = np.sqrt(X**2 + Y**2)
            PHI = np.arctan2(Y, X)
            bc = np.linspace(i, pointsall.shape[0]-1, pointsall.shape[0], dtype=int)
        else:
            try:
                nodes, elementi, bc, mesher = MKE_solver().read(mreza)
                elementi = elementi[:,1:]-1
                vozlisca = np.empty((elementi.shape[0],3,2), dtype=object)
                for i in range(len(vozlisca)):
                    vozlisca[i] = nodes[elementi[i],1:]
                bc = bc['boundary_nodes']-1 #index!
                n_node,X,Y = np.transpose(nodes)
                pointsall = np.transpose([X,Y])
                R = np.sqrt(X**2 + Y**2)
                PHI = np.arctan2(Y, X)
                generate_mesh = False
            except: raise ValueError('Mreža ni podana pravilno!')
        figsize=(2*Nc*2+0.5, Nr*2)
        # unique_pts, counts = np.unique(pointsall, axis=0, return_counts=True)
        # duplicate_indices = np.where(counts > 1)[0]
        # duplicate_points = unique_pts[duplicate_indices]
        # print("Duplicate points:")
        # print(duplicate_points)            
    else: raise ValueError('Opna ni podana pravilno!')

    if generate_mesh: # triangulacija mreže
        triang = Delaunay(pointsall)
        elementi = triang.simplices #indeksi vozlišč po KE, shape(št.el, št.vozlišč na KE=3)
        vozlisca = pointsall[triang.simplices] #koordinate vozlišč po KE, shape(št.el, št.vozlišč na KE=3, dimenzije=2)
    N_el = elementi.shape[0] #št. elementov
    N_voz = pointsall.shape[0] #št. vozlišč

    # zlaganje
    try:
        A=np.zeros((N_voz,N_voz)) #leva stran - togostna matrika
        B=np.zeros((N_voz,N_voz)) #desna stran - masna matrika
    except: #if N_voz == 64592: #memory error
        A=np.zeros((N_voz,N_voz), dtype=np.float32) #leva stran - togostna matrika
        B=np.zeros((N_voz,N_voz), dtype=np.float32) #desna stran - masna matrika
    def d(a,b): return (a[0]-b[0])**2 + (a[1]-b[1])**2
    for i in range(N_el):
        voz_index = elementi[i]
        voz = vozlisca[i]
        (x1,y1),(x2,y2),(x3,y3) = voz
        S = 0.5 * np.abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))
        K = np.asarray([
            [d(voz[1], voz[2]), (x2-x3)*(x3-x1)+(y2-y3)*(y3-y1), (x3-x2)*(x2-x1)+(y3-y2)*(y2-y1)],
            [(x1-x3)*(x3-x2)+(y1-y3)*(y3-y2), d(voz[0], voz[2]), (x3-x1)*(x1-x2)+(y3-y1)*(y1-y2)],
            [(x1-x2)*(x2-x3)+(y1-y2)*(y2-y3), (x2-x1)*(x1-x3)+(y2-y1)*(y1-y3), d(voz[1], voz[0])]])
        A[np.transpose([voz_index]), voz_index] += K/4/S
        B[np.transpose([voz_index]), voz_index] += S/12*(np.eye(3)+1)
    # zero_rows = np.where(np.all(B == 0, axis=1))[0]
    # print("Zero rows in B:", zero_rows)

    # RP
    A[bc,:] = 0
    A[:,bc] = 0
    # B[bc,:] = 0
    # B[:,bc] = 0
    A[bc,bc] = 1e6 #povečana togost -> pomik = 0, ne prevelika vrednost - stabilnost!

    # lastne rešitve
    try: val, vec = eigsh(A, min([Nr*Nc*2,N_voz]), B, sigma=0) #sparse hitreje!
    except: val, vec = eigh(A,B)
    val = np.real(val)
    vec = np.real(vec)
    k = np.sqrt(val) #lastne vrednosti
    arg = np.argsort(k)
    k = k[arg]
    vec = vec[:,arg]
    if '.dat' in mreza and 'mesh_1.dat' not in mreza and 'mesh_01.dat' not in mreza: #čudno težava le ansys gostejše mreže ???
        k*=1e3
        vec*=1e3
    t2=timer()
    print('\nšt. vozlišč:',N_voz)

    if plot:
        print('št. vozlišč:',N_voz)
        print('št. elementov:',N_el)
        print('t = ',t2-t1, ' s')

        plt.triplot(pointsall[:,0], pointsall[:,1], elementi, c='gray')
        plt.plot(pointsall[:,0], pointsall[:,1], 'k.')
        plt.plot(pointsall[bc,0], pointsall[bc,1], 'rx')
        plt.axis('equal')
        plt.xlabel('$x$')
        plt.ylabel('$y$')
        plt.tight_layout()
        # plt.savefig(f'./slike/mreza_{opna}_{mreza}_N{N_voz}.pdf')
        plt.show()

        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(Nr,Nc)
        ax = gs.subplots(sharex=True, sharey=True)
        vmax,vmin=np.max(vec[:,:Nr*Nc]), np.min(vec[:,:Nr*Nc])
        i=j=0
        m,n=1,1
        k_analit=[]
        k_napaka = []
        v_napaka = []
        for l in range(Nr*Nc):
            if opna == 'polkrog':
                _,_,u,ki = analit(m,n, r=R, phi=PHI)
                k_analit.append(ki)
                k_napaka.append(np.abs(k[l]-ki))
                v_napaka.append(np.linalg.norm(vec[:,l]-u))
            elif opna == 'kvadrat':
                # u = np.sin(k[l]*X.ravel())*np.sin(k[l]*Y.ravel())
                _,_,u,ki = analit_sq(m,n, X,Y)
                k_analit.append(ki)
                k_napaka.append(np.abs(k[l]-ki))
                v_napaka.append(np.linalg.norm(vec[:,l]-u.ravel()))
            cs=ax[j,i].tricontourf(X.ravel(),Y.ravel(),vec[:,l], vmin=vmin, vmax=vmax, levels=10, cmap='turbo')
            # cs=ax[j,i].tricontourf(X.ravel(),Y.ravel(),np.abs(vec[:,l])-np.abs(u), levels=10, cmap='turbo')
            if opna == 'f': ax[j,i].plot(X.ravel()[bc], Y.ravel()[bc], 'k.')
            ax[j,i].grid(False)
            ax[j,i].axis('off')
            ax[j,i].set_xticks([])
            ax[j,i].set_yticks([])
            ax[j,i].set_title(f'$k_{{{i+1}}}$ = {k[l]:.3g}', fontsize=11)
            i+=1
            n+=1
            if n+1>m: n=1; m+=1
            if i==Nc: i=0; j+=1
        cbar_ax = fig.add_axes([0.86, 0, 0.07, 1])
        cbar_ax.axis('off')
        cbar=fig.colorbar(cs, ax=cbar_ax, orientation='vertical')
        cbar.ax.tick_params(labelsize=11)
        cbar.set_label('u', fontsize=12)
        plt.savefig(f'./slike/last_obl_{opna}_{mreza}_N{N_voz}.pdf')
        plt.show()

        plt.plot(k[:30],'kx-', label='MKE')
        plt.plot(k_analit[:30],'r.', label='analitično')
        # if opna == 'polkrog': plt.plot(k_g[:30],'g.', label='Galerkin')
        plt.legend()
        plt.xlabel('i')
        plt.ylabel('lastna vrednost $k$')
        plt.minorticks_on()
        plt.tight_layout()
        # plt.savefig(f'./slike/last_vred_{opna}_{mreza}_N{N_voz}.pdf')
        plt.show()

        plt.plot(k_napaka, 'kx-', label='$|k_i-k_{ex}|$')
        plt.plot(v_napaka, 'rx-', label='norma razlik lastne oblike')
        plt.yscale('log')
        plt.legend()
        plt.xlabel('index')
        plt.ylabel('absolutna napaka')
        plt.minorticks_on()
        plt.tight_layout()
        # plt.savefig(f'./slike/napaka_{opna}_{mreza}_N{N_voz}.pdf')
        plt.show()

        # m=n=1
        # for i in range(2*N):
        #     _,_,u,ki = analit(m,n, r=R, phi=PHI) #le polkrog!
        #     print('i = ',i,', u_max = ',np.max(np.abs(vec[:,i])))
        #     # tocke=plt.tricontourf(X.ravel(),Y.ravel(),vec[:,i],levels=10) #delaunay triangulacija!
        #     # tocke=plt.tricontourf(Xg.ravel(),Yg.ravel(),vec_g[:,i],levels=10) #delaunay triangulacija!
        #     tocke=plt.tricontourf(X.ravel(),Y.ravel(),np.abs(vec[:,i])-np.abs(u),levels=10) #delaunay triangulacija!
        #     # tocke=plt.scatter(points[:,0],points[:,1], c=vec[:,i], cmap='turbo', linewidths=0)
        #     plt.axis('equal')
        #     plt.grid(False)
        #     plt.xlabel('$x$')
        #     plt.ylabel('$y$')
        #     plt.title(f'$k_{{{i+1}}}$ = {k[i]:.3g}')
        #     plt.colorbar(tocke, label='u')
        #     plt.show()
        #     n+=1
        #     if n+1>m: n=1; m+=1
    return k, t2-t1, N_voz

def richardson1(A, r=2):
    """definicija po RDTM le za prve (zadnje boljše?) TRI vrednosti!, 
    abolutno - ni težav z oscilacijami, vendar ne upoštevamo predznaka konvergence!"""
    # p = np.log(np.abs((A[0]-A[1]) / (A[1]-A[2]))) / np.log(r) #prvi
    # A_ex = A[2] + (A[1]-A[2]) / (r**p - 1)
    p = np.log(np.abs((A[-1]-A[-2]) / (A[-2]-A[-3]))) / np.log(r) #zadnji
    A_ex = A[-2] + (A[-2]-A[-3]) / (r**p - 1)
    return A_ex,p
def richardson2(A_k, r=2, p=2):
    '''za skalarne veličine, red zgostitve r=2, predpostavljen red konvergence p=2'''
    A_k = np.array(A_k)
    for k in range(1,len(A_k)):
        A_k_old = np.copy(A_k)
        A_k = np.zeros(A_k_old.shape[0]-1)
        for i in range(len(A_k_old)-1):
            data_1 = A_k_old[i] #data_2: 2x gostejša mreža od data_1
            data_2 = A_k_old[i+1]
            A_k[i] = 1/(r**(p+k)-1)*(r**(p+k)*data_2 - data_1)
    return A_k[0]
def linfit(x, y):
    '''linearna regresija'''
    A = np.vstack([x, np.ones(len(x))]).T
    k, n = np.linalg.lstsq(A, y, rcond=None)[0]
    return k,n

N = 30
opna = 'polkrog' # 'polkrog' ali 'kvadrat' ali 'f'
mreza = 'cartesian' #za polkrog ('polar', 'random polar', 'cartesian', 'random cartesian')
Nc,Nr=5,4 #št. stolpcev in vrstic za prikaz lastnih rešitev

# solver(N, opna, mreza, Nc, Nr, plot=True)
# solver(N=30, opna='kvadrat', mreza='', Nc=Nc, Nr=Nr, plot=True)
# solver(N=100, opna='f', mreza='', Nc=Nc, Nr=Nr, plot=True)
# solver(N, 'polkrog', './mesh/mesh_05.inp', Nc, Nr, plot=True)
# solver(N, 'polkrog', 'mesh/mesh_80mm.dat', Nc, Nr, plot=True)

Nc,Nr=2,2 #št. stolpcev in vrstic za prikaz lastnih rešitev
calc(M=5, N=5, eigs=Nc*Nr, to_plot_u=True, Nc=Nc,Nr=Nr)
input('konec')

# f opna
#################################################################################################
N_list = np.linspace(10, 50, 5, dtype=int)
N_list = np.linspace(10, 130, 10, dtype=int)
# N_list = np.array([2**i for i in range(4,7)]) #3 za richardson
t, k, N_voz = [], [], []
N_eig = 7 #le prvih nekaj lastne vrednosti!

print(' ')
for i,N in enumerate(N_list): #Nx=Ny=N
    print('lastna mreža N = ',N, ', ',i+1,'/',len(N_list), end='\r')
    kn,tn,Nv = solver(N, 'polkrog', 'cartesian', Nc, Nr)
    k.append(kn[:N_eig]) #le prve tri lastne vrednosti!
    t.append(tn)
    N_voz.append(Nv)
N_voz = np.asarray(N_voz)
k = np.asarray(k)
t = np.asarray(t)

fig,(ax1) = plt.subplots(1,1, figsize=(7,4), sharex=True)
c=['r','g','b','c','m','orange','y']
p=[]
for i in range(k.shape[1]):
    k1,n1 = linfit(np.log10(1/N_voz), np.log10(np.abs(k[:,0])))
    A1,p1 = richardson1(k[:,i])
    Ak1 = richardson2(k[:,i], p=p1)
    p.append(p1)
    
    ax1.semilogx(N_voz, k[:,i], 'x-',c=c[i], label=f'{i+1}.')
    ax1.axhline(A1, color=c[i], linestyle='--')
    ax1.axhline(Ak1, color=c[i], linestyle='-.')
ax1.plot(Ak1, color='k', linestyle='--', label='Richardson 1')
ax1.plot(Ak1, color='k', linestyle='-.', label='Richardson 2')
ax1.minorticks_on()
ax1.legend(title='$k_i$', loc='upper left', bbox_to_anchor=(1.05, 1))
ax1.set_xlabel('$N$')
ax1.set_ylabel('$k$')
plt.tight_layout()
# plt.savefig(f'./slike/last_vred.pdf')
plt.show()

p= np.array(p)
fig,ax1 = plt.subplots(1,1, figsize=(6,4))
ax1.bar(np.arange(1,k.shape[1]+0.1,dtype=int), np.abs(p), color='gray')
ax1.grid(False)
ax1.set_xlabel('$k_i$')
ax1.set_ylabel('red konvergence $|p|$')
plt.tight_layout()
# plt.savefig(f'./slike/konvergence.pdf')
plt.show()

fig,ax1 = plt.subplots(1,1, figsize=(6,4), sharex=True)
ax1.plot(N_voz, t, 'kx-')
ax1.minorticks_on()
ax1.set_xlabel('$N$')
ax1.set_ylabel('$t$ [s]')
plt.tight_layout()
# plt.savefig(f'./slike/time.pdf')
plt.show()

# polkrog
########################################################################
k={'cartesian':[], 'random cartesian':[], 'polar':[], 'random polar':[], 'abaqus':[], 'ansys':[], 'galerkin':[]}
t={'cartesian':[], 'random cartesian':[], 'polar':[], 'random polar':[], 'abaqus':[], 'ansys':[], 'galerkin':[]}
N_voz={'cartesian':[], 'random cartesian':[], 'polar':[], 'random polar':[], 'abaqus':[], 'ansys':[], 'galerkin':[]}

def name1(i): return f'/home/ravni/Desktop/Modelska_II/mesh/mesh_{i}.inp'
def name2(i): return f'/home/ravni/Desktop/Modelska_II/mesh/mesh_{i}mm.dat'

# N_list = np.linspace(10, 35, 10, dtype=int)
N_list = np.linspace(10, 130, 10, dtype=int)
#zgoščevanje zadnja memory error !!
# abaq_data = [name1(i) for i in ['1', '05', '025', '0125', '01', '005']]
# ansys_data = ['/home/ravni/Desktop/Modelska_II/mesh/mesh_1.dat', name2('80'), name2('40'), name2('20'), '/home/ravni/Desktop/Modelska_II/mesh/mesh_01.dat', name2('5')]
abaq_data = [name1(i) for i in ['1', '05', '025', '0125', '01']] #zgoščevanje
ansys_data = ['/home/ravni/Desktop/Modelska_II/mesh/mesh_1.dat', name2('80'), name2('40'), name2('20'), '/home/ravni/Desktop/Modelska_II/mesh/mesh_01.dat'] #zgoščevanje

N_list = np.array([2**i for i in range(4,7)]) #3 za richardson
abaq_data = [name1(i) for i in ['1', '05', '025']]
ansys_data = [name2('80'), name2('40'), name2('20')]

print(' ')
for i,N in enumerate(N_list): #Nx=Ny=N
    print('lastna mreža N = ',N, ', ',i+1,'/',len(N_list), end='\r')
    kn,tn,Nv = solver(N, 'polkrog', 'cartesian', Nc, Nr)
    k['cartesian'].append(kn[:3])
    t['cartesian'].append(tn)
    N_voz['cartesian'].append(Nv)
    kn,tn,Nv = solver(N, 'polkrog', 'random cartesian', Nc, Nr)
    k['random cartesian'].append(kn[:3])
    t['random cartesian'].append(tn)
    N_voz['random cartesian'].append(Nv)
    kn,tn,Nv = solver(N, 'polkrog', 'polar', Nc, Nr)
    k['polar'].append(kn[:3])
    t['polar'].append(tn)
    N_voz['polar'].append(Nv)
    kn,tn,Nv = solver(N, 'polkrog', 'random polar', Nc, Nr)
    k['random polar'].append(kn[:3])
    t['random polar'].append(tn)
    N_voz['random polar'].append(Nv)
# galerkin
print('')
N_max = 5
for N in range(2,N_max+1,1):
    for M in range(1,N):
        print('Galerkin M = ',M, ', N = ',N ,'/',N_max, end='\r')
        try:
            kn,_,tn = calc(M=5, N=5, eigs=3, to_calc_u=True)
            k['galerkin'].append(kn[:3])
            t['galerkin'].append(tn)
            N_voz['galerkin'].append(M*N)
        except: 
            input('Galerkin napaka pri M=',M,', N=',N, end='\r')
            pass
print('')
for i in range(len(abaq_data)):
    print('Abaqus',i+1,'/', len(abaq_data), end='\r')
    kn,tn,Nv = solver(N, 'polkrog', abaq_data[i], Nc, Nr)
    k['abaqus'].append(kn[:3])
    t['abaqus'].append(tn)
    N_voz['abaqus'].append(Nv)
print('')
for i in range(len(ansys_data)): #čudno skaliranje z 1000 !!??
    print('Ansys',i+1,'/', len(ansys_data), end='\r')
    kn,tn,Nv = solver(N, 'polkrog', ansys_data[i], Nc, Nr)
    k['ansys'].append(kn[:3])
    t['ansys'].append(tn)
    N_voz['ansys'].append(Nv)
print('')
for key in k: k[key] = np.array(k[key])
for key in N_voz: N_voz[key] = np.array(N_voz[key])
for key in t: t[key] = np.array(t[key])
k_ex = [analit(m=1,n=1)[3], analit(m=2,n=1)[3], analit(m=3,n=1)[3]] #analitične lastne vrednosti

# plotting
########################################################################
A_richardson = {'cartesian':[], 'random cartesian':[], 'polar':[], 'random polar':[], 'abaqus':[], 'ansys':[], 'galerkin':[]}
conv_order = {'cartesian':[], 'random cartesian':[], 'polar':[], 'random polar':[], 'abaqus':[], 'ansys':[], 'galerkin':[]} #p
kn_fit = {'cartesian':[], 'random cartesian':[], 'polar':[], 'random polar':[], 'abaqus':[], 'ansys':[], 'galerkin':[]}
fig,(ax1,ax2,ax3) = plt.subplots(3,1, figsize=(8, 10), sharex=True)
c=['r','g','b','c','m','orange','y']
for i, key in enumerate(k.keys()):
    A1,p1 = richardson1(k[key][:,0]) #le zgostitev x2 !!!
    A2,p2 = richardson1(k[key][:,1])
    A3,p3 = richardson1(k[key][:,2])
    Ak1 = richardson2(k[key][:,0], p=p1)
    Ak2 = richardson2(k[key][:,1], p=p2)
    Ak3 = richardson2(k[key][:,2], p=p3)
    A_richardson[key].append([[A1,A2,A3], [Ak1,Ak2,Ak3]])
    conv_order[key].append([p1,p2,p3])

    k1,n1 = linfit(np.log10(1/N_voz[key]), np.log10(np.abs(k[key][:,0]-k_ex[0])))
    k2,n2 = linfit(np.log10(1/N_voz[key]), np.log10(np.abs(k[key][:,1]-k_ex[1])))
    k3,n3 = linfit(np.log10(1/N_voz[key]), np.log10(np.abs(k[key][:,2]-k_ex[2])))
    kn_fit[key].append([[k1,k2,k3], [n1,n2,n3]])
    try:
        ax1.axhline(k_ex[0], color='k', linestyle='--')
        ax1.semilogx(N_voz[key], k[key][:,0], 'x-', color=c[i], label=key)
        ax2.axhline(k_ex[1], color='k', linestyle='--')
        ax2.semilogx(N_voz[key], k[key][:,1], 'x-', color=c[i])
        ax3.axhline(k_ex[2], color='k', linestyle='--')
        ax3.semilogx(N_voz[key], k[key][:,2], 'x-', color=c[i])
    except: pass
ax1.plot(N_list[0], k['cartesian'][0,0], 'k--', label='analitično') #dummy
ax1.minorticks_on()
ax1.legend()
ax3.set_xlabel('$N$')
ax1.set_ylabel('$k_1$')
ax2.set_ylabel('$k_2$')
ax3.set_ylabel('$k_3$')
plt.tight_layout()
plt.savefig(f'./slike/last_vred.pdf')
plt.show()

for key in k: A_richardson[key] = np.array(A_richardson[key])
for key in k: conv_order[key] = np.array(conv_order[key])
for key in k: kn_fit[key] = np.array(kn_fit[key])

fig,(ax1,ax2) = plt.subplots(2,1, figsize=(8, 20), sharex=True)
bars=np.array([0,1,2])
bar_width=0.1
for i, key in enumerate(k.keys()):
    ax1.bar(bars, np.abs(kn_fit[key][0][0]), label=key, width=bar_width, color=c[i], edgecolor='k')
    ax2.bar(bars, np.abs(kn_fit[key][0][1]), label=key, width=bar_width, color=c[i], edgecolor='k')
    bars=bars+bar_width
ax1.minorticks_on()
ax2.minorticks_on()
ax1.legend()
ax1.grid(False, axis='x')
ax2.grid(False, axis='x')
ax2.set_xticks([0.25,1.25,2.25], labels=['$k_1$', '$k_2$', '$k_3$'])
ax2.set_xlabel('lastne vresnosti')
ax1.set_ylabel('$k$')
ax2.set_ylabel('$n$')
plt.tight_layout()
plt.savefig(f'./slike/linfit.pdf')
plt.show()

fig,(ax1) = plt.subplots(1,1, figsize=(8, 10), sharex=True)
bars=np.array([0,1,2])
bar_width=0.1
for i, key in enumerate(k.keys()):
    ax1.bar(bars, np.abs(conv_order[key][0]), label=key, width=bar_width, color=c[i], edgecolor='k')
    bars=bars+bar_width
ax1.minorticks_on()
ax1.legend()
ax1.grid(False, axis='x')
ax1.set_xticks([0.25,1.25,2.25], labels=['$k_1$', '$k_2$', '$k_3$'])
ax1.set_xlabel('lastne vresnosti')
ax1.set_ylabel('red konvergence $p$')
plt.tight_layout()
plt.savefig(f'./slike/conv_order.pdf')
plt.show()

fig,(ax1,ax2) = plt.subplots(2,1, figsize=(8, 20), sharex=True)
bars=np.array([0,1,2])
bar_width=0.1
for i, key in enumerate(k.keys()):
    ax1.bar(bars, np.abs(A_richardson[key][0][0]-k_ex)/k_ex, label=key, width=bar_width, color=c[i], edgecolor='k')
    ax2.bar(bars, np.abs(A_richardson[key][0][1]-k_ex)/k_ex, label=key, width=bar_width, color=c[i], edgecolor='k')
    bars=bars+bar_width
ax1.minorticks_on()
ax2.minorticks_on()
ax1.legend()
ax1.set_yscale('log')
ax2.set_yscale('log')
ax1.grid(False, axis='x')
ax2.grid(False, axis='x')
ax2.set_xticks([0.25,1.25,2.25], labels=['$k_1$', '$k_2$', '$k_3$'])
ax2.set_xlabel('lastne vresnosti')
ax1.set_ylabel('Richardson1 relativni pogrešek')
ax2.set_ylabel('Richardson2 relativni pogrešek')
plt.tight_layout()
plt.savefig(f'./slike/richardson_err.pdf')
plt.show()

fig,(ax1,ax2,ax3) = plt.subplots(3,1, figsize=(8, 10), sharex=True)
for i, key in enumerate(k.keys()):
    try:
        ax1.loglog(1/N_voz[key], np.abs(k[key][:,0] - k_ex[0]), 'x-', color=c[i], label=key)
        ax2.loglog(1/N_voz[key], np.abs(k[key][:,1] - k_ex[1]), 'x-', color=c[i])
        ax3.loglog(1/N_voz[key], np.abs(k[key][:,2] - k_ex[2]), 'x-', color=c[i])

        nx=np.array([1/N_voz[key][0]*1.5, 1/N_voz[key][-1]*0.5])
        ax1.loglog(nx, 10**(kn_fit[key][0][0][0]*np.log10(nx)+kn_fit[key][0][1][0]), '--', color=c[i])
        ax2.loglog(nx, 10**(kn_fit[key][0][0][1]*np.log10(nx)+kn_fit[key][0][1][1]), '--', color=c[i])
        ax3.loglog(nx, 10**(kn_fit[key][0][0][2]*np.log10(nx)+kn_fit[key][0][1][2]), '--', color=c[i])
    except: pass
ax1.minorticks_on()
ax1.legend()
ax3.set_xlabel('$h^*$')
ax1.set_ylabel('$|k_1-k_{ex}|$')
ax2.set_ylabel('$|k_2-k_{ex}|$')
ax3.set_ylabel('$|k_3-k_{ex}|$')
plt.tight_layout()
plt.savefig(f'./slike/last_vred_err_h.pdf')
plt.show()

fig,ax1 = plt.subplots(1,1, figsize=(8, 10), sharex=True)
for i, key in enumerate(k.keys()):
    try: ax1.plot(N_voz[key], t[key], 'x-', color=c[i], label=key)
    except: pass
ax1.minorticks_on()
ax1.legend()
ax1.set_xlabel('$N$')
ax1.set_ylabel('$t$ [s]')
plt.tight_layout()
plt.savefig(f'./slike/time.pdf')
plt.show()
