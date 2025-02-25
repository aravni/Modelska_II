# -*- coding: utf-8 -*-
# spektralno reševanje poissonove enačbe na kvadratni domeni

# %%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import scipy.fft as fft
from timeit import default_timer as timer

plt.style.use('aljaz')
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

def gauss_seidel(A,q,N=1,tol=1e-6):
    D = np.diag(A)
    u_n = np.zeros_like(q)
    n = 0
    residual = tol+1

    #iteracije
    while n < N and residual > tol:
        n += 1
        # u_n = np.zeros_like(u)
        # Dirichlet RP
        # for i in range(Nx): #spodaj
        #     u_n[i] = q[i]
        # for j in range(1,Ny): #levo
        #     u_n[j*Nx] = q[j*Nx]
        # for i in range(Nx): #gor
        #     u_n[-Nx+i] = q[-Nx+i]
        # for j in range(Ny): #desno
        #     u_n[j*Nx+Nx-1] = q[j*Nx+Nx-1]

        for j in range(1,Ny-1):
            for i in range(1,Nx-1):
                ij_index = j*Nx + i
                u_n[ij_index] = 1/D[ij_index] * (q[ij_index] + 1/dx**2*(u_n[ij_index+1] + u_n[ij_index-1])
                                                             + 1/dy**2*(u_n[ij_index+Nx] + u_n[ij_index-Nx]))
        residual = np.linalg.norm(A@u_n-q)
        # u = u_n
    if n == N: print('Max iterations reached!')
    return u_n, residual, n

def SOR(A,q,omega=1,N=1,tol=1e-6):
    D = np.diag(A)
    u = np.zeros_like(q)
    n = 0
    residual = tol+1

    #iteracije
    while n < N and residual > tol:
        n += 1
        # u_n = np.zeros_like(u)
        # # Dirichlet RP -> se prepiše od prejšnje iteracije

        for j in range(1,Ny-1):
            for i in range(1,Nx-1):
                ij_index = j*Nx + i
                u_gs = 1/D[ij_index] * (q[ij_index] + 1/dx**2*(u[ij_index+1] + u[ij_index-1])
                                                             + 1/dy**2*(u[ij_index+Nx] + u[ij_index-Nx]))
                u[ij_index] = u[ij_index] + omega*(u_gs - u[ij_index])
        #residual = np.linalg.norm(u_n-u)
        residual = np.linalg.norm(A@u-q) #lahko počasi!
        # u = u_n #ni treba shranjevat dveh iteraij
    if n == N: print('Max iterations reached!')
    return u, residual, n

a=b=1 #širina=višina
N=51 #št. računskih točk po osi
h=a/(N-1)
ratio=100 #ratio=1 - homogeno

def g_homog(x,y,r):
    '''desna stran PDE homogena opna'''
    return -np.ones_like(r) #homogena opna
def g_kvad(x,y,r, ratio=2):
    '''desna stran PDE kvadratno območje odebelitve'''
    l = 0.5*a #stranica kvadrata
    origin = (a/2,b/2) #izhodišče
    origin = (3*a/4, 3*b/4) #desno zgoraj
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

x = np.linspace(0,a,N)
y = np.linspace(0,b,N)
r=np.zeros((len(x),len(y)))
X,Y = np.meshgrid(x,y)

#diskretne vrednosti desne strani v matrični obliki
g_jk = g_homog(x,y,r)
ime= 'homogena'

# g_jk = g_kvad(x,y,r,ratio=ratio)
# ime= 'nehomogena-kvadrat na sredini'
# ime= 'nehomogena-kvadrat na robu'

# g_jk = g_f(x,y,r,ratio=ratio)
# ime= 'f-opna'

plt.figure(figsize=(4,4))
plt.imshow(g_jk,origin='lower', cmap='binary_r')
# plt.scatter(X,Y,c=g_jk.flatten())
plt.axis('equal')
# plt.axis('off')
plt.xticks([])
plt.yticks([])
plt.xlim(0,N-2)
plt.ylim(0,N-2)
# plt.title('oblika opne: '+ime)
plt.savefig(f'./slike/{ime}')
plt.show()

u_homog = lambda u_max: u_max*16*(X-X**2)*(Y-Y**2) #?analitična? rešitev za homogeno kvadratno domeno
# k = np.arange(0,N*N,1,dtype=int)
# u_homog = g_jk.flatten()/np.sinh(k*np.pi)*np.sin(k*np.pi*X.flatten())*np.sinh(k*np.pi*Y.flatten()) #NI OK
# print(1-np.max(u_homog(1)))
# u_homog = u_homog(1)#.reshape((N,N))
# sc=plt.contourf(x,y,u_homog,cmap='turbo')
# cbar=plt.colorbar(sc,label='u')
# plt.axis('equal')
# plt.show()

# %%    1. 2x 1D sinusna FFT
t = timer()
G_mk = fft.dst(g_jk, axis=0, type=1) #pol-transformiran prostor
G_mn = fft.dst(G_mk, axis=1, type=1) #transformiran prostor

m=n=np.arange(1,1+N,1,dtype=int)
m,n=np.meshgrid(m,n)

U_mn = h*h/2 * G_mn / (np.cos(np.pi*m/N)+np.cos(np.pi*n/N)-2) #vektorsko
U_jn = fft.idst(U_mn, axis=0, type=1) #1D sinusna FT
u_jk = fft.idst(U_jn, axis=1, type=1) #nazaj v fiz prostor
t_fft = timer()-t
u_1d_sin = u_jk
u_h = u_homog(np.max(u_jk))
print('1. 1D FFT.',t_fft,',',np.linalg.norm(u_jk-u_h), np.max(u_1d_sin))
title = '1D sin fft'

# %%    2. 2D sinusna FFT
t = timer()
G_mn = fft.dstn(g_jk, type=1) #transformiran prostor

m=n=np.arange(1,1+N,1,dtype=int)
m,n=np.meshgrid(m,n)

U_mn = h*h/2 * G_mn / (np.cos(np.pi*m/N)+np.cos(np.pi*n/N)-2) #vektorsko
u_jk = fft.idstn(U_mn, type=1) #nazaj v fiz prostor
t_fft = timer()-t
u_2d_sin = u_jk
u_h = u_homog(np.max(u_jk))
print('2. 2D FFT.',t_fft,',',np.linalg.norm(u_jk-u_h), np.max(u_2d_sin)) #očitno sistem enačb(A@u-g) ni dobro merilo!!!!!!
title = '2D sin fft'

# U=np.zeros((N,N)) #vsote (v trasformiranem prostoru) -> počasno!!
# for m1 in range(1,N+1):
#     for n1 in range(1,N+1):
#         U[m1-1,n1-1] = h**2/2 * G_mn[m1-1,n1-1] / (np.cos(np.pi*m1/N)+np.cos(np.pi*n1/N)-2)
# u_naiv= fft.idstn(U) #nazaj v fiz prostor
# print(np.sum(np.abs(u_naiv-u_jk))) #isto

# %%    3. zvezna FT, naivno seštevanje vsot (v realni domeni), POČASNO
t = timer()
m=n=np.arange(1,1+N,1,dtype=int)
m,n=np.meshgrid(m,n)
A_mn = -4*(1-(-1)**m)*(1-(-1)**n)/(np.pi**2*m*n*((m*np.pi/a)**2+(n*np.pi/b)**2)) #le za homogeno: g(x,y) = 1 !!!
u = np.zeros((N,N))
for m1 in range(1,N+1):
    for n1 in range(1,N+1):
        u -= A_mn[m1-1,n1-1] * np.sin(m1*np.pi*X/a) * np.sin(n1*np.pi*Y/b)
u_jk = u
t_ft = timer()-t
u_zvezno = u_jk
u_h = u_homog(np.max(u_jk))
print('3. zvezna FT:',t_ft,'s ,',np.linalg.norm(u_jk-u_h), np.max(u_zvezno))
title = '2D zvezna ft'

# %%    4. kombinirana metoda
t = timer()
pas_sr = -(4-2*np.cos(m*np.pi/N))
pas_sp = pas_zg = np.ones((pas_sr.size-1))

G_m_k = fft.dst(g_jk, axis=1, type=1)
desna = h**2*G_m_k
# RP: povsod na robu Dirichlet NIČ! -> že v DST!
Trp=0 #Dirichlet RP (spl)
print(desna.shape,pas_sp.shape,pas_sr.shape, pas_zg.shape)

U_m_k = TDMAsolver(pas_sp, pas_sr, pas_zg, desna)

u_jk = fft.idst(U_m_k, axis=1, type=1)
u_jk += Trp
t_komb = timer()-t
u_komb = u_jk
u_h = u_homog(np.max(u_jk))
print('4. kombinirana:',t_komb,'s ,',np.linalg.norm(u_jk-u_h), np.max(u_komb))
title = 'kombinirana'

# plt.contourf(X,Y,u_jk)
# plt.colorbar()
# plt.show()

# %%    5. SOR
N_iter = 500 #max iter
tol = 1e-3 #toleranca rešitve (leva - desna stran), 10^-6 je prehudo!
omega = 2/(1+np.pi/N) #ocena
t = timer()
title = 'SOR'

# matrični zapis sistema enačb
Nx=Ny=N #konst. mreža
dx=dy=h
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

u_jk, residual_sor, n_iter_sor = SOR(A,g, omega, N=N_iter, tol=tol)
u_jk = u_jk.reshape((N,N))
t_sor = timer()-t
u_sor = u_jk
u_h = u_homog(np.max(u_jk))
print('5. SOR:',t_sor,'s ,', np.linalg.norm(u_jk-u_h), f'N_iter = {n_iter_sor}, residual = {residual_sor}, {np.max(u_sor)}')

# %%    6. Gauss-Seidel -> zelo počasi!
N_iter = 10*N_iter #max iter - 10x več od SOR!
# tol = 1e-3 #toleranca rešitve (leva - desna stran) - isto kot SOR
t = timer()
title = 'Gauss-Seidel'

u_jk, residual_gs, n_iter_gs = gauss_seidel(A,g, N=N_iter, tol=tol)
u_jk = u_jk.reshape((N,N))
t_gs = timer()-t
u_gs = u_jk
u_h = u_homog(np.max(u_jk))
print('6. Gauss-Seidel:',t_gs,'s ,', np.linalg.norm(u_jk-u_h), f'N_iter = {n_iter_gs}, residual = {residual_gs}, {np.max(u_gs)}')

t=timer()
u_lu = np.linalg.solve(A,g).reshape((N,N))
u_h = u_homog(np.max(u_lu))
print(f'7. Eksaktno: {timer()-t} s, {np.linalg.norm(u_lu-u_h)}, {np.max(u_lu)}')

# %%    plot
pointsize=3
u=u_jk.reshape((N,N))
max_u = np.max(u)

# sc=plt.contourf(x,y,u_jk,cmap='turbo')
# #sc=plt.scatter(x,y,c=u,s=pointsize)
# cbar=plt.colorbar(sc,label='u')
# cbar_ticks = list(cbar.get_ticks()[:-1])  # Get existing ticks
# cbar_ticks.append(max_u)  # Add max(u) to the ticks
# cbar.set_ticks(cbar_ticks)
# plt.axis('equal')
# plt.xlabel('$x$')
# plt.ylabel('$y$')
# plt.title(title+f', $N={N}$')
# plt.tight_layout()
# #plt.savefig('./u_'+ime+'.pdf')
# plt.show()

#vse
fig, ax = plt.subplots(1,7, figsize=(20,4), gridspec_kw={'wspace':0.1})
vmin = min(u_2d_sin.min(), u_komb.min(), u_sor.min(), u_gs.min(), u_zvezno.min())
vmax = max(u_2d_sin.max(), u_komb.max(), u_sor.max(), u_gs.max(), u_zvezno.max())
v_max= max([vmax,abs(vmin)])
vmin /= v_max
vmax /= v_max
levels = np.round(np.linspace(vmin, vmax, 11),1)

ax_i = ax[0]
cont=ax_i.contourf(X, Y, u_2d_sin/v_max, cmap='turbo', vmin=vmin, vmax=vmax, levels=levels)
ax_i.axis('equal')
ax_i.set_title(f'2D sin FFT\n', wrap=True)
ax_i.axis('off')

ax_i = ax[1]
cont=ax_i.contourf(X, Y, u_komb/v_max, cmap='turbo', vmin=vmin, vmax=vmax, levels=levels)
ax_i.axis('equal')
ax_i.set_title(f'kombinirana\n', wrap=True)
ax_i.axis('off')

ax_i = ax[2]
cont=ax_i.contourf(X, Y, u_sor/v_max, cmap='turbo', vmin=vmin, vmax=vmax, levels=levels)
ax_i.axis('equal')
ax_i.set_title(f'SOR\n$N_{{iter}}={n_iter_sor}$', wrap=True)
ax_i.axis('off')

ax_i = ax[3]
cont=ax_i.contourf(X, Y, u_gs/v_max, cmap='turbo', vmin=vmin, vmax=vmax, levels=levels)
ax_i.axis('equal')
ax_i.set_title(f'Gauss-Seidel\n$N_{{iter}}={n_iter_gs}$', wrap=True)
ax_i.axis('off')

ax_i = ax[4]
cont=ax_i.contourf(X, Y, u_lu/v_max, cmap='turbo', vmin=vmin, vmax=vmax, levels=levels)
ax_i.set_title(f'eksaktna rešitev\nsistema enačb', wrap=True)
ax_i.axis('equal')
ax_i.axis('off')

# OK le za homogena opna!
ax_i = ax[5]
ax_i.contourf(X, Y, u_zvezno/v_max, cmap='turbo', vmin=vmin, vmax=vmax, levels=levels)
ax_i.set_title(f'zvezna FT\n')
ax_i.axis('equal')
ax_i.axis('off')

ax_i = ax[6]
cont=ax_i.contourf(X, Y, u_1d_sin/v_max, cmap='turbo', vmin=vmin, vmax=vmax, levels=levels)
ax_i.axis('equal')
ax_i.set_title(f'1D sin FFT\n', wrap=True)
ax_i.axis('off')

cbar_ax = fig.add_axes([0.88, 0.175, 0.05, 0.7])
cbar_ax.axis('off')
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
sm = cm.ScalarMappable(cmap='turbo', norm=norm)
sm.set_array([])
cbar=fig.colorbar(sm, ax=cbar_ax)
cbar.set_label('$\\frac{u}{u_{max}}$', fontsize=20)
cbar.set_ticks(levels)

fig.suptitle(f'št. računskih točk $N={N*N}$')
plt.subplots_adjust(top=0.75)
plt.tight_layout(rect=[0, 0, 0.92, 0.9])
plt.savefig(f'./slike/metode_{ime}_r{ratio}', bbox_inches='tight') #shranjevanje OK, prikaz ne!
plt.show()

# %% max poves
fig=plt.figure(figsize=(10,4))
gs=fig.add_gridspec(1,2, wspace=0)
ax = gs.subplots(sharey=True)
ax[0].plot(x, u_1d_sin[N//2,:], label='1D sin FFT')
ax[0].plot(x, u_2d_sin[N//2,:], '--',label='2D sin FFT')
ax[0].plot(x, u_zvezno[N//2,:], label='zvezna FT')
ax[0].plot(x, u_komb[N//2,:], label='kombinirana')
ax[0].plot(x, u_sor[N//2,:], label='SOR')
ax[0].plot(x, u_gs[N//2,:], label='Gauss-Seidel')
ax[0].plot(x, u_lu[N//2,:], label='eksaktno')
ax[0].set_xlabel('$x$')
ax[0].set_ylabel('$u$')
ax[0].minorticks_on()

ax[1].plot(x, u_1d_sin[:,N//2], label='1D sin FFT')
ax[1].plot(x, u_2d_sin[:,N//2], '--', label='2D sin FFT')
ax[1].plot(x, u_zvezno[:,N//2], label='zvezna FT')
ax[1].plot(x, u_komb[:,N//2], label='kombinirana')
ax[1].plot(x, u_sor[:,N//2], label='SOR')
ax[1].plot(x, u_gs[:,N//2], label='Gauss-Seidel')
ax[1].plot(x, u_lu[:,N//2], label='eksaktno')
ax[1].set_xlabel('$y$')
ax[1].minorticks_on()
ax[1].legend(loc='upper left', bbox_to_anchor=(1,1))
plt.tight_layout()
plt.savefig(f'./slike/prerez_{ime}_r{ratio}')
plt.show()

u_max = np.array([np.max(u_1d_sin), np.max(u_2d_sin), np.max(u_zvezno), np.max(u_komb), np.max(u_sor), np.max(u_gs), np.max(u_lu)])
plt.figure(figsize=(8,4))
plt.plot(u_max, 'kx')
plt.xticks(np.arange(7),['1D sin\nFFT', '2D sin\nFFT', 'zvezna\nFT', 'kombinirana', 'SOR', 'Gauss-Seidel', 'eksaktno'])
plt.ylabel('$u_{max}$')
plt.minorticks_on()
plt.tight_layout()
plt.savefig(f'./slike/Umax_{ime}_r{ratio}')
plt.show()