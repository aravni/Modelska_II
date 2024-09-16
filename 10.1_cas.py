# -*- coding: utf-8 -*-
# spektralno reševanje poissonove enačbe na kvadratni domeni

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

def SOR(A,q,N,omega=1,N_iter=1,tol=1e0):
    D = np.diag(A)
    u = np.zeros_like(q)
    n = 0
    residual = tol+1
    Nx=Ny=N #konst. mreža
    dx=dy=a/N
    #iteracije
    while n < N_iter and residual > tol:
        n += 1
        for j in range(1,Ny-1):
            for i in range(1,Nx-1):
                ij_index = j*Nx + i
                u_gs = 1/D[ij_index] * (q[ij_index] + 1/dx**2*(u[ij_index+1] + u[ij_index-1])
                                                             + 1/dy**2*(u[ij_index+Nx] + u[ij_index-Nx]))
                u[ij_index] = u[ij_index] + omega*(u_gs - u[ij_index])
        residual = np.linalg.norm(A@u-q)
    if n == N_iter: print(f'Solution did not converged to {tol:.1e}')
    return u, residual, n

def g_homog(x,y,r):
    '''desna stran PDE homogena opna'''
    return -np.ones_like(r) #homogena opna
def g_kvad(x,y,r, ratio=2):
    '''desna stran PDE kvadratno območje odebelitve'''
    l = 0.5*a #stranica kvadrata
    origin = (a/2,b/2) #izhodišče
    origin = (3*a/4, 3*b/4)
    for i in range(r.shape[0]):
        for j in range(r.shape[1]):
            if abs(x[i]-origin[0]) < l/2 and abs(y[j]-origin[1]) < l/2:
                r[j,i] = ratio
            else: r[j,i] = 1
    return -r #nehomogena - debelejša v kvadratu
def g_f(x,y,r, ratio=2):
    '''desna stran PDE oblika f'''
    for i in range(r.shape[0]):
        for j in range(r.shape[1]):
            if a/4 < x[i] < a/2:
                r[j,i] = ratio
            elif y[j] > 3*b/4 and x[i] > a/4:
                r[j,i] = ratio
            elif b/4 < y[j] < b/2 and x[i] < 3*a/4:
                r[j,i] = ratio
            else: r[j,i] = 1
    return -r

u_homog = lambda u_max: u_max*16*(X-X**2)*(Y-Y**2) #?analitična? rešitev za homogeno kvadratno domeno

def sin_fft(N):
    t = timer()
    U_mn = h*h/2 * G_mn / (np.cos(np.pi*m/N)+np.cos(np.pi*n/N)-2) #vektorsko
    u_jk = fft.idstn(U_mn) #nazaj v fiz prostor
    t_fft = timer()-t
    u_h = u_homog(np.max(u_jk))
    residual = np.linalg.norm(u_jk-u_h)
    return u_jk, t_fft, residual

def zvezna_ft(N):
    t = timer()
    A_mn = -4*(1-(-1)**m)*(1-(-1)**n)/(np.pi**2*m*n*((m*np.pi/a)**2+(n*np.pi/b)**2)) #le za homogeno: g(x,y) = 1 !!!
    u = np.zeros((N,N))
    for m1 in range(1,N+1):
        for n1 in range(1,N+1):
            u -= A_mn[m1-1,n1-1] * np.sin(m1*np.pi*X/a) * np.sin(n1*np.pi*Y/b)
    u_jk = u
    t_ft = timer()-t
    u_h = u_homog(np.max(u_jk))
    residual = np.linalg.norm(u_jk-u_h)
    return u_jk, t_ft, residual

def kombinirana(N):
    t = timer()
    pas_sr = -(4-2*np.cos(m*np.pi/N))
    pas_sp = pas_zg = np.ones((pas_sr.size-1))

    G_m_k = fft.dst(g_jk,axis=1)
    desna = h**2*G_m_k

    U_m_k = TDMAsolver(pas_sp, pas_sr, pas_zg, desna)
    u_jk = fft.idst(U_m_k,axis=1)
    t_komb = timer()-t
    u_h = u_homog(np.max(u_jk))
    residual = np.linalg.norm(u_jk-u_h)
    return u_jk, t_komb, residual

def SOR_solver(N):
    N_iter = 1000 #max iteracij
    omega = 2/(1+np.pi/N) #ocena
    t = timer()

    # matrični zapis sistema enačb
    Nx=Ny=N #konst. mreža
    dx=dy=a/N #param a se ne spreminja!
    A = np.zeros((Nx*Ny, Nx*Ny))
    g = -g_jk.flatten()
    # DE
    for j in range(1,Ny-1):
        for i in range(1,Nx-1):
            ij_index = j*(Nx) + i
            ij = 2*(1/dx**2 + 1/dy**2)
            ij_dx = -1/dx**2 #- 1/(2*i*dx**2) #zadnji člen = polarni dodatek!
            ij_dy = -1/dy**2
            #center
            A[ij_index, ij_index] += ij
            # okolica
            A[ij_index, ij_index-1] += ij_dx #levo
            A[ij_index, ij_index+1] += ij_dx #desno
            A[ij_index, ij_index-Nx] += ij_dy #dol
            A[ij_index, ij_index+Nx] += ij_dy #gor
    # RP Dirichlet
    rp = 0
    for i in range(Nx):
        #spodaj
        A[i, i] = 1
        g[i] = rp
        #zgoraj
        A[-Nx+i, -Nx+i] = 1
        g[-Nx+i] = rp
    for j in range(Ny):
        #levo
        A[j*Nx, j*Nx] = 1
        g[j*Nx] = rp
        #desno
        A[j*Nx+Nx-1, j*Nx+Nx-1] = 1
        g[j*Nx+Nx-1] = rp

    #u_jk, residual, n_iter = gauss_seidel(A,g, N=N_iter) #zelo počasi!
    u_jk, residual, n_iter = SOR(A,g, N, omega, N_iter)
    u_jk = u_jk.reshape((N,N))
    t_sor = timer()-t
    u_h = u_homog(np.max(u_jk))
    residual = np.linalg.norm(u_jk-u_h)
    return u_jk, t_sor, residual, n_iter

######################################################################
a=b=1. #širina=višina
ratio = 10 #razmerje mas - hehomog
#test
# N=10
# h=a/N
# x = np.linspace(0,1,N)
# y = np.linspace(0,1,N)
# r=np.zeros((len(x),len(y)))
# X,Y = np.meshgrid(x,y)
# g_jk = g_f(x,y,r,ratio=ratio)
# G_mn = fft.dstn(g_jk)
# m=n=np.arange(1,1+N,1,dtype=int)
# m,n=np.meshgrid(m,n)
# u= sin_fft(N)
# plot=plt.contourf(x,y,u[0],cmap='turbo')
# plt.colorbar(plot,label='u')
# plt.axis('equal')
# plt.show()

N_list_1=np.logspace(1,2.2,30,dtype=int) #(1,2.2,10) #za SOR in zvezno
N_list=np.logspace(1,4,60,dtype=int) #(1,4,50)
u = {'sin fft':[], 'zvezna ft':[], 'kombinirana':[], 'sor':[]}
čas = {'sin fft':[], 'zvezna ft':[], 'kombinirana':[], 'sor':[]}
residuali = {'sin fft':[], 'zvezna ft':[], 'kombinirana':[], 'sor':[]}
n_iter = []

i=1
for N in N_list:
    print(f'{i}/{len(N_list)}',end='\r')
    h=a/N
    i+=1

    x = np.linspace(0,a,N)
    y = np.linspace(0,b,N)
    r=np.zeros((len(x),len(y)))
    X,Y = np.meshgrid(x,y)

    #diskretne vrednosti desne strani v matrični obliki
    # g_jk = g_homog(x,y,r)
    # g_jk = g_kvad(x,y,r,ratio=ratio)
    g_jk = g_f(x,y,r,ratio=ratio)

    #transformiran prostor
    G_mn = fft.dstn(g_jk)
    m=n=np.arange(1,1+N,1,dtype=int)
    m,n=np.meshgrid(m,n)

    u_r, t_r, r_r = sin_fft(N)
    u['sin fft'].append(np.max(u_r))
    čas['sin fft'].append(t_r)
    residuali['sin fft'].append(r_r)

    u_r, t_r, r_r = kombinirana(N)
    u['kombinirana'].append(np.max(u_r))
    čas['kombinirana'].append(t_r)
    residuali['kombinirana'].append(r_r)

print('\n')
i=1
for N in N_list_1:
    print(f'{i}/{len(N_list_1)}',end='\r')
    h=a/N
    i+=1

    x = np.linspace(0,a,N)
    y = np.linspace(0,b,N)
    r=np.zeros((len(x),len(y)))
    X,Y = np.meshgrid(x,y)

    #diskretne vrednosti desne strani v matrični obliki
    # g_jk = g_homog(x,y,r)
    # g_jk = g_kvad(x,y,r,ratio=ratio)
    g_jk = g_f(x,y,r,ratio=ratio)

    #transformiran prostor
    G_mn = fft.dstn(g_jk)
    m=n=np.arange(1,1+len(G_mn),1,dtype=int)
    m,n=np.meshgrid(m,n)

    # zvezna le za homogeno opno!!!
    # u_r, t_r, r_r = zvezna_ft(N)
    # u['zvezna ft'].append(np.max(u_r))
    # čas['zvezna ft'].append(t_r)
    # residuali['zvezna ft'].append(r_r)

    u_r, t_r, r_r, n_r = SOR_solver(N)
    u['sor'].append(np.max(u_r))
    čas['sor'].append(t_r)
    residuali['sor'].append(r_r)
    n_iter.append(n_r)

coef = np.polyfit(np.log(N_list), np.log(čas['sin fft']), deg=1)
O_fft = coef[0]
lin_fft = np.polyval(coef,np.log(N_list))
# coef = np.polyfit(np.log(N_list_1), np.log(čas['zvezna ft']), deg=1)
# O_ft = coef[0]
# lin_ft = np.polyval(coef,np.log(N_list_1))
coef = np.polyfit(np.log(N_list), np.log(čas['kombinirana']), deg=1)
O_komb = coef[0]
lin_komb = np.polyval(coef,np.log(N_list))
coef = np.polyfit(np.log(N_list_1), np.log(čas['sor']), deg=1)
O_sor = coef[0]
lin_sor = np.polyval(coef,np.log(N_list_1))

fig,ax=plt.subplots()
ax.plot(N_list, čas['sin fft'],'C0.',label=f'2D sin fft $\mathcal{{O}}$($N^{{{O_fft:.2f}}}$)')
ax.plot(N_list, čas['kombinirana'],'C1.',label=f'kombinirana $\mathcal{{O}}$($N^{{{O_komb:.2f}}}$)')
# ax.plot(N_list_1, čas['zvezna ft'],'C2.',label=f'zvezna ft $\mathcal{{O}}$($N^{{{O_ft:.2f}}}$)')
ax.plot(N_list_1, čas['sor'],'C3.',label=f'SOR $\mathcal{{O}}$($N^{{{O_sor:.2f}}}$)')

ax1=ax.twinx()
ax1.plot(N_list, lin_fft, 'C0-')
ax1.plot(N_list, lin_komb, 'C1-')
# ax1.plot(N_list_1, lin_ft, 'C2-')
ax1.plot(N_list_1, lin_sor, 'C3-')
ax1.axis('off')

ax.grid(which='minor', c='gray',alpha=0.25)
ax.grid(which='major',c='k')
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend()
ax.set_xlabel('$N$')
ax.set_ylabel('$t$ [s]')
plt.tight_layout()
plt.savefig('./slike/cas_N_opna.pdf')
plt.show()

plt.plot(N_list_1, n_iter,'.-')
plt.grid(which='minor', c='gray',alpha=0.25)
plt.grid(which='major',c='k')
plt.xlabel('$N$')
plt.ylabel('$N_{iter}$')
plt.xscale('log')
plt.title('število iteracij SOR za doseganje pogreška 1')
plt.savefig('./slike/Niter_N_opna.pdf')
plt.show()

plt.plot(N_list, residuali['sin fft'],'C0.-',label=f'2D sin fft')
plt.plot(N_list, residuali['kombinirana'],'C1.-',label=f'kombinirana ')
# plt.plot(N_list_1, residuali['zvezna ft'],'C2.-',label=f'zvezna ft ')
plt.plot(N_list_1, residuali['sor'],'C3.-',label=f'SOR')
plt.grid(which='minor', c='gray',alpha=0.25)
plt.grid(which='major',c='k')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.xlabel('$N$')
plt.ylabel('residual')
plt.title('odstopanje od polinomske rešitve')
plt.savefig('./slike/residuali_N_kvad_opna.pdf')
plt.show()

plt.plot(N_list, u['sin fft'],'C0.-',label=f'2D sin fft')
plt.plot(N_list, u['kombinirana'],'C1.-',label=f'kombinirana ')
# plt.plot(N_list_1, u['zvezna ft'],'C2.-',label=f'zvezna ft ')
plt.plot(N_list_1, u['sor'],'C3.-',label=f'SOR')
plt.grid(which='minor', c='gray',alpha=0.25)
plt.grid(which='major',c='k')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.xlabel('$N$')
plt.ylabel('$u_{max}$')
plt.title('konvergenca max povesa')
plt.tight_layout()
plt.savefig('./slike/konvergenca_kvad_opna.pdf')
plt.show()