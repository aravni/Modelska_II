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

def solver(N, opna, mreza, Nc, Nr, plot=False):
    t1=timer()
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
        else: raise ValueError('Mreža ni podana pravilno!')
        figsize=(2*Nc*2+0.5, Nr*2)
        # unique_pts, counts = np.unique(pointsall, axis=0, return_counts=True)
        # duplicate_indices = np.where(counts > 1)[0]
        # duplicate_points = unique_pts[duplicate_indices]
        # print("Duplicate points:")
        # print(duplicate_points)

        if mreza == 'Galerkin': #TODO
            N_r, N_phi = 10,10
            dr_g=1/(N_r-1)
            d_phi_g=2*np.pi/(N_phi-1)
            r_g=np.linspace(dr_g,1,N_r)
            phi_g=np.linspace(0,np.pi,N_phi)
            R_g,PHI_g=np.meshgrid(r_g,phi_g)
            R_g=np.append(R_g,0) #le 1x izhodišče!
            PHI_g=np.append(PHI_g,0)
            Xg=R_g*np.cos(PHI_g)
            Yg=R_g*np.sin(PHI_g)

            Ag = np.zeros((N_r*N_phi,N_r*N_phi))
            Bg = np.zeros((N_r*N_phi,N_r*N_phi))
            for m in range(N_phi):
                A_m = np.zeros((N_r,N_r))
                B_m = np.zeros((N_r,N_r))
                m+=1
                for k in range(1,N_r+1):
                    for l in range(1,N_r+1):
                        A_m[k-1,l-1] = np.pi/2 * (k*l/(2*m+k+l) - (2*k*l+k+l)/(2*m+k+l+1) + (k+1)*(l+1)/(2*m+k+l+2))
                        B_m[k-1,l-1] = np.pi/2 * (1/(2*m+k+l+2) - 2/(2*m+k+l+3) + 1/(2*m+k+l+4))
                m-=1
                Ag[(m)*N_r:(m+1)*N_r, (m)*N_r:(m+1)*N_r] = A_m
                Bg[(m)*N_r:(m+1)*N_r, (m)*N_r:(m+1)*N_r] = B_m
            # lastne rešitve
            val_g, vec_g = eig(Ag,Bg)
            val_g = np.real(val_g)
            vec_g = np.real(vec_g)
            k_g = np.sort(np.sqrt(val_g)) #lastne vrednosti
            t2=timer()
            print('Galerkin: t = ',t2-t1, ' s') 
    else: raise ValueError('Opna ni podana pravilno!')

    # triangulacija mreže
    triang = Delaunay(pointsall)
    elementi = triang.simplices #indeksi vozlišč po KE, shape(št.el, št.vozlišč na KE=3)
    vozlisca = pointsall[triang.simplices] #koordinate vozlišč po KE, shape(št.el, št.vozlišč na KE=3, dimenzije=2)

    N_el = elementi.shape[0] #št. elementov
    N_voz = pointsall.shape[0] #št. vozlišč

    # zlaganje
    A=np.zeros((N_voz,N_voz)) #leva stran - togostna matrika
    B=np.zeros((N_voz,N_voz)) #desna stran - masna matrika
    def d(a,b): return (a[0]-b[0])**2 + (a[1]-b[1])**2
    def phi_approx(m,n): return r**(m+n) * (1-r) * np.sin(m*phi)
    for i in range(N_el):
        voz = vozlisca[i]
        voz_index = elementi[i]
        (x1,y1),(x2,y2),(x3,y3) = voz
        S = 0.5 * np.abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))
        K = np.array([
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
    t2=timer()

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
                u = np.sin(k[l]*X.ravel())*np.sin(k[l]*Y.ravel())
                k_napaka.append(None)
                v_napaka.append(np.linalg.norm(vec[:,l]-u))
            cs=ax[j,i].tricontourf(X.ravel(),Y.ravel(),vec[:,l], vmin=vmin, vmax=vmax, levels=10, cmap='turbo')
            # cs=ax[j,i].tricontourf(X.ravel(),Y.ravel(),np.abs(vec[:,l])-np.abs(u), levels=10, cmap='turbo')
            if opna == 'f': ax[j,i].plot(X.ravel()[bc], Y.ravel()[bc], 'k.')
            ax[j,i].grid(False)
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
        plt.show()

        plt.plot(k[:30],'kx-', label='MKE')
        if opna == 'polkrog': plt.plot(k_analit[:30],'r.', label='analitično')
        if opna == 'polkrog': plt.plot(k_g[:30],'g.', label='Galerkin')
        plt.legend()
        plt.xlabel('i')
        plt.ylabel('lastna vrednost $k$')
        plt.minorticks_on()
        plt.show()

        plt.plot(k_napaka, 'kx-', label='$|k_i-k_{ex}|$')
        plt.plot(v_napaka, 'rx-', label='norma razlik lastne oblike')
        plt.yscale('log')
        plt.legend()
        plt.xlabel('index')
        plt.ylabel('absolutna napaka')
        plt.minorticks_on()
        plt.show()

        m=n=1
        for i in range(2*N):
            _,_,u,ki = analit(m,n, r=R, phi=PHI) #le polkrog!
            print('i = ',i,', u_max = ',np.max(np.abs(vec[:,i])))
            # tocke=plt.tricontourf(X.ravel(),Y.ravel(),vec[:,i],levels=10) #delaunay triangulacija!
            # tocke=plt.tricontourf(Xg.ravel(),Yg.ravel(),vec_g[:,i],levels=10) #delaunay triangulacija!
            tocke=plt.tricontourf(X.ravel(),Y.ravel(),np.abs(vec[:,i])-np.abs(u),levels=10) #delaunay triangulacija!
            # tocke=plt.scatter(points[:,0],points[:,1], c=vec[:,i], cmap='turbo', linewidths=0)
            plt.axis('equal')
            plt.grid(False)
            plt.xlabel('$x$')
            plt.ylabel('$y$')
            plt.title(f'$k_{{{i+1}}}$ = {k[i]:.3g}')
            plt.colorbar(tocke, label='u')
            plt.show()
            n+=1
            if n+1>m: n=1; m+=1
    return k, t2-t1

N = 30
opna = 'polkrog' # 'polkrog' ali 'kvadrat' ali 'f'
mreza = 'cartesian' #za polkrog ('polar', 'random polar', 'cartesian', 'random cartesian')
Nc,Nr=5,4 #št. stolpcev in vrstic za prikaz lastnih rešitev
# solver(N, opna, mreza, Nc, Nr, plot=True)

# polkrog
########################################################################
k={'cartesian':[], 'random cartesian':[], 'polar':[], 'random polar':[]}
t={'cartesian':[], 'random cartesian':[], 'polar':[], 'random polar':[]}
N_list = np.linspace(10,30,10, dtype=int) #[2**2, 2**3, 2**4, 2**5]
for N in N_list:
    print('N = ',N)
    kn,tn = solver(N, 'polkrog', 'cartesian', Nc, Nr)
    k['cartesian'].append(kn[:3])
    t['cartesian'].append(tn)
    kn,tn = solver(N, 'polkrog', 'random cartesian', Nc, Nr)
    k['random cartesian'].append(kn[:3])
    t['random cartesian'].append(tn)
    kn,tn = solver(N, 'polkrog', 'polar', Nc, Nr)
    k['polar'].append(kn[:3])
    t['polar'].append(tn)
    kn,tn = solver(N, 'polkrog', 'random polar', Nc, Nr)
    k['random polar'].append(kn[:3])
    t['random polar'].append(tn)
for key in k: k[key] = np.array(k[key])
k_ex = [analit(m=1,n=1)[3], analit(m=2,n=1)[3], analit(m=3,n=1)[3]] #analitične lastne vrednosti
# TODO richardson exstrapolacija, red konvergence!!! - RDT

fig,(ax1,ax2,ax3) = plt.subplots(3,1, figsize=(8, 10), sharex=True)
c=['r','g','b','c','m','y']
for i, key in enumerate(k.keys()):
    ax1.axhline(k_ex[0], color='k', linestyle='--')
    ax1.plot(N_list, k[key][:,0], 'x-', color=c[i], label=key)
    ax2.axhline(k_ex[1], color='k', linestyle='--')
    ax2.plot(N_list, k[key][:,1], 'x-', color=c[i])
    ax3.axhline(k_ex[2], color='k', linestyle='--')
    ax3.plot(N_list, k[key][:,2], 'x-', color=c[i])
ax1.plot(N_list[0], k['cartesian'][0,0], 'k--', label='analitično') #dummy
ax1.minorticks_on()
ax1.legend()
ax3.set_xlabel('$N$')
ax1.set_ylabel('$k_1$')
ax2.set_ylabel('$k_2$')
ax3.set_ylabel('$k_3$')
plt.tight_layout()
plt.show()

fig,(ax1,ax2,ax3) = plt.subplots(3,1, figsize=(8, 10), sharex=True)
c=['r','g','b','c','m','y']
for i, key in enumerate(k.keys()):
    ax1.loglog(N_list, np.abs(k[key][:,0] - k_ex[0]), 'x-', color=c[i], label=key)
    ax2.loglog(N_list, np.abs(k[key][:,1] - k_ex[1]), 'x-', color=c[i])
    ax3.loglog(N_list, np.abs(k[key][:,2] - k_ex[2]), 'x-', color=c[i])
ax1.minorticks_on()
ax1.legend()
ax3.set_xlabel('$N$')
ax1.set_ylabel('$k_1$')
ax2.set_ylabel('$k_2$')
ax3.set_ylabel('$k_3$')
plt.tight_layout()
plt.show()

fig,ax1 = plt.subplots(1,1, figsize=(8, 10), sharex=True)
c=['r','g','b','c','m','y']
for i, key in enumerate(k.keys()):
    ax1.plot(N_list, t[key], 'x-', color=c[i], label=key)
ax1.minorticks_on()
ax1.legend()
ax1.set_xlabel('$N$')
ax1.set_ylabel('$t$ [s]')
plt.tight_layout()
plt.show()