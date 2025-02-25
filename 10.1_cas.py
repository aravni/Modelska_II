# -*- coding: utf-8 -*-
# spektralno reševanje poissonove enačbe na kvadratni domeni
# z različnimi porazdelitvami mase

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import scipy.fft as fft
from timeit import default_timer as timer
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from tqdm import tqdm
import os

plt.style.use('aljaz')
def pause():
    input("Press the <ENTER> key to continue...")

def parallel(data):
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

    def sin_fft_2d(N):
        t = timer()
        G_mn = fft.dstn(g_jk)
        U_mn = h*h/2 * G_mn / (np.cos(np.pi*m/N)+np.cos(np.pi*n/N)-2) #vektorsko
        u_jk = fft.idstn(U_mn) #nazaj v fiz prostor
        t_fft = timer()-t
        u_h = u_homog(np.max(u_jk))
        residual = np.linalg.norm(u_jk-u_h)
        return u_jk, t_fft, residual
    
    def sin_fft_1d(N):
        t = timer()
        G_mk = fft.dst(g_jk, axis=0) #pol-transformiran prostor
        G_mn = fft.dst(G_mk, axis=1) #transformiran prostor
        m=n=np.arange(1,1+N,1,dtype=int)
        m,n=np.meshgrid(m,n)
        U_mn = h*h/2 * G_mn / (np.cos(np.pi*m/N)+np.cos(np.pi*n/N)-2) #vektorsko
        U_jn = fft.idst(U_mn, axis=0) #1D sinusna FT
        u_jk = fft.idst(U_jn, axis=1) #nazaj v fiz prostor
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

    N,a,ratio = data
    b=a #širina=višina -> kvadrat
    h=a/N

    x = np.linspace(0,a,N)
    y = np.linspace(0,b,N)
    r=np.zeros((len(x),len(y)))
    X,Y = np.meshgrid(x,y)

    #diskretne vrednosti desne strani v matrični obliki
    g_jk = g_homog(x,y,r)
    # g_jk = g_kvad(x,y,r,ratio=ratio)
    # g_jk = g_f(x,y,r,ratio=ratio)

    #transformiran prostor
    m=n=np.arange(1,1+N,1,dtype=int)
    m,n=np.meshgrid(m,n)

    u_r1, t_r1, r_r1 = sin_fft_2d(N)
    u_r2, t_r2, r_r2 = kombinirana(N)
    u_r4, t_r4, n_r4 = sin_fft_1d(N)
    if N < 10**2.2: #max N za SOR in zevezno -> počasni!
        u_r3, t_r3, r_r3, n_r3 = SOR_solver(N)  
        u_r5, t_r5, n_r5 = zvezna_ft(N)
    else:
        u_r3, t_r3, r_r3, n_r3 = [None],None,None,None  #TODO vpliv None na sledeče izračune!!
        u_r5, t_r5, n_r5 = [None],None,None
    return u_r1, t_r1, r_r1, u_r2, t_r2, r_r2, u_r3, t_r3, r_r3, n_r3, u_r4, t_r4, n_r4, u_r5, t_r5, n_r5

######################################################################
if __name__ == '__main__':
    nthrds=8
    a=b=1. #širina=višina
    ratio = 10 #razmerje mas - nehomog
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

    os.environ["OMP_NUM_THREADS"] = str(1)  # OpenMP
    os.environ["MKL_NUM_THREADS"] = str(1)  # Intel MKL
    os.environ["OPENBLAS_NUM_THREADS"] = str(1)  # OpenBLAS
    os.environ["NUMEXPR_NUM_THREADS"] = str(1)  # NumExpr

    # N_list_1=np.logspace(1,2.2,30,dtype=int) #(1,2.2,10) #za SOR in zvezno
    # N_list=np.logspace(1,4,60,dtype=int) #(1,4,50)
    u = {'sin fft':[], '1d sin':[], 'zvezna':[], 'kombinirana':[], 'sor':[]}
    čas = {'sin fft':[], '1d sin':[], 'zvezna':[], 'kombinirana':[], 'sor':[]}
    residuali = {'sin fft':[], '1d sin':[], 'zvezna':[], 'kombinirana':[], 'sor':[]}
    n_iter = []

    N_list=np.logspace(1,2.3,10,dtype=int) #test

    # Function to update the progress bar
    def update_progress(*args):
        progress_bar.update()
    # Initialize the progress bar
    progress_bar = tqdm(total=len(N_list), desc='Processing')
    progress_bar.update()

    tasks = [(N,a,ratio) for N in N_list]
    with ProcessPoolExecutor(nthrds) as executor:
        # Create a list of futures
        futures = [executor.submit(parallel, task) for task in tasks]
        for future in futures:
            future.add_done_callback(update_progress)
        # Wait for all futures to complete and collect results
        results = [future.result() for future in futures]
        # Close the progress bar
        progress_bar.close()
    
    for result in results:
        u['sin fft'].append(np.max(result[0]))
        čas['sin fft'].append(result[1])
        residuali['sin fft'].append(result[2])
        u['kombinirana'].append(np.max(result[3]))
        čas['kombinirana'].append(result[4])
        residuali['kombinirana'].append(result[5])
        u['sor'].append(np.max(result[6]))
        čas['sor'].append(result[7])
        residuali['sor'].append(result[8])
        n_iter.append(result[9])
        u['1d sin'].append(np.max(result[10]))
        čas['1d sin'].append(result[11])
        residuali['1d sin'].append(result[12])
        u['zvezna'].append(np.max(result[13]))
        čas['zvezna'].append(result[14])
        residuali['zvezna'].append(result[15])
    input('end')

    coef = np.polyfit(np.log(N_list), np.log(čas['sin fft']), deg=1)
    O_fft = coef[0]
    lin_fft = np.polyval(coef,np.log(N_list))
    coef = np.polyfit(np.log(N_list), np.log(čas['1d sin']), deg=1)
    O_fft_1d = coef[0]
    lin_fft_1d = np.polyval(coef,np.log(N_list))
    coef = np.polyfit(np.log(N_list), np.log(čas['zvezna']), deg=1)
    O_ft = coef[0]
    lin_ft = np.polyval(coef,np.log(N_list))
    coef = np.polyfit(np.log(N_list), np.log(čas['kombinirana']), deg=1)
    O_komb = coef[0]
    lin_komb = np.polyval(coef,np.log(N_list))
    coef = np.polyfit(np.log(N_list), np.log(čas['sor']), deg=1)
    O_sor = coef[0]
    lin_sor = np.polyval(coef,np.log(N_list))

    fig,ax=plt.subplots()
    ax.plot(N_list, čas['sin fft'],'C0.',label=f'2D sin FFT $\mathcal{{O}}$($N^{{{O_fft:.2f}}}$)')
    ax.plot(N_list, čas['kombinirana'],'C1.',label=f'kombinirana $\mathcal{{O}}$($N^{{{O_komb:.2f}}}$)')
    ax.plot(N_list, čas['zvezna'],'C2.',label=f'zvezna FT $\mathcal{{O}}$($N^{{{O_ft:.2f}}}$)')
    ax.plot(N_list, čas['sor'],'C3.',label=f'SOR $\mathcal{{O}}$($N^{{{O_sor:.2f}}}$)')
    ax.plot(N_list, čas['1d sin'],'C4.',label=f'1D sin FFT $\mathcal{{O}}$($N^{{{O_fft_1d:.2f}}}$)')

    ax1=ax.twinx()
    ax1.plot(N_list, lin_fft, 'C0-')
    ax1.plot(N_list, lin_komb, 'C1-')
    ax1.plot(N_list, lin_ft, 'C2-')
    ax1.plot(N_list, lin_sor, 'C3-')
    ax1.plot(N_list, lin_fft_1d, 'C4-')
    ax1.axis('off')

    plt.minorticks_on()
    ax.grid(which='both')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    ax.set_xlabel('$N$')
    ax.set_ylabel('$t$ [s]')
    plt.tight_layout()
    # plt.savefig('./slike/cas_homog.pdf')
    plt.show()

    # N_iter SOR
    plt.plot(N_list, n_iter,'k.-')
    plt.grid(which='minor', c='gray',alpha=0.25)
    plt.grid(which='major',c='k')
    plt.xlabel('$N$')
    plt.ylabel('$N_{iter}$')
    plt.xscale('log')
    plt.title('število iteracij SOR za doseganje pogreška 1')
    # plt.savefig('./slike/Niter_homog.pdf')
    plt.show()

    # residual glede na homogeno polinomsko rešitev
    plt.plot(N_list, residuali['sin fft'],'C0.-',label=f'2D sin FFT')
    plt.plot(N_list, residuali['kombinirana'],'C1.-',label=f'kombinirana ')
    plt.plot(N_list, residuali['zvezna'],'C2.-',label=f'zvezna FT ')
    plt.plot(N_list, residuali['sor'],'C3.-',label=f'SOR')
    plt.plot(N_list, residuali['1d sin'],'C4.--',label=f'1D sin FFT')
    plt.xscale('log')
    # plt.yscale('log')
    plt.minorticks_on()
    ax.grid(which='both')
    plt.legend()
    plt.xlabel('$N$')
    plt.ylabel('residual')
    plt.title('odstopanje od polinomske rešitve')
    # plt.savefig('./slike/residuali_homog.pdf')
    plt.show()

    # max poves
    plt.plot(N_list, u['sin fft'],'C0.-',label=f'2D sin FFT')
    plt.plot(N_list, u['kombinirana'],'C1.-',label=f'kombinirana ')
    plt.plot(N_list, u['zvezna'],'C2.-',label=f'zvezna FT ')
    plt.plot(N_list, u['sor'],'C3.-',label=f'SOR')
    plt.plot(N_list, u['1d sin'],'C4.--',label=f'1D sin FFT')

    plt.xscale('log')
    # plt.yscale('log')
    plt.minorticks_on()
    ax.grid(which='both')
    plt.legend()
    plt.xlabel('$N$')
    plt.ylabel('$u_{max}$')
    plt.title('konvergenca max povesa')
    plt.tight_layout()
    # plt.savefig('./slike/konvergenca_Umax_homog.pdf')
    plt.show()