# -*- coding: utf-8 -*-
# 4. Hartree-Fock metoda
# 1. del: ionizacijska energija He atoma

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.integrate import solve_ivp
from scipy import optimize as opt
from scipy.integrate import simps,trapz, cumtrapz, cumulative_simpson
from timeit import default_timer as timer
from numba import njit, jit
from matplotlib.colors import ListedColormap, BoundaryNorm
from alignYaxes import *
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator

plt.style.use('aljaz')

@njit
def k_fun(x,epsilon,phi,z): #z = št. elektronov
    return 2*z/x + 2*phi + epsilon
def DE(x,y,epsilon,phi,z):
    R, z = y.reshape((2,len(y)//2))
    dRdx = z
    dzdx = k_fun(x,epsilon,phi,z) * R
    return np.append(dRdx,dzdx)
def numerov(x, y0, e, phi, z):
    '''reševanje DE oblike [d²/dx² + k²(x)] y(x) = 0,
    funkcija k(x) naj bo že kvadrirana. Pazi na časovni korak!'''
    h = x[1]-x[0]
    y = np.zeros_like(x)

    y[0] = y0[0] #ZP
    y[1] = y0[0]+h*y0[1] #ZP aprox iz naklona
    # y[1] = y0[1] #ZP 2. točka

    for i in range(1,len(x)-1):
        a_2 = 1 - h**2*5/12*k_fun(x[i],e,phi[i],z)
        a_1 = 1 + h**2*1/12*k_fun(x[i-1],e,phi[i-1],z)
        a_3 = 1 + h**2*1/12*k_fun(x[i+1],e,phi[i+1],z)
        y[i+1] = (2*a_2*y[i] - a_1*y[i-1]) / a_3
    return y
def runge_kutta_4(f, t, y0, e, phi, z):
    """
    Metoda Runge-Kutta 4. reda za reševanje diferencialne enačbe: y' = f(t, y)
    
    :param f: funkcija, ki jo kličemo s parametroma t in y in vrne 
              vrednost prvega odvoda
    :param t: ekvidistanten časovni vektor oz. neodvisna spremenljivka
    :param y0: začetna vrednost
    :return y: funkcijske vrednosti.
    """   
    y = np.zeros(shape=(len(t),len(y0)))
    y[0] = y0
    h = t[1]-t[0]

    i=0
    for ti in t[1:]:
        k0 = h*f(ti, y[i], e, phi[i],z)
        k1 = h*f(ti + h/2.0, y[i] + k0/2.0, e, phi[i],z)
        k2 = h*f(ti + h/2.0, y[i] + k1/2.0, e, phi[i],z)
        k3 = h*f(ti + h, y[i] + k2, e, phi[i],z)
        y[i+1] = y[i] + (k0 + 2.0*k1 + 2.0*k2 + k3)/6.0
        i+=1
    return np.transpose(y)
def calc_potential(R,x):
    """Calculate the Hartree-Fock potential."""
    phi = np.zeros_like(x)
    for i in range(1,len(x)):
        # phi[i] = -1/x[i] * trapz(R[:i]**2, x[:i]) - trapz(R[i:]**2/x[i:], x[i:])
        phi[i] = -1/x[i] * simps(R[:i]**2, x[:i]) - simps(R[i:]**2/x[i:], x[i:])
    # phi[1:] = -1/x[1:] * cumtrapz(R**2,x) - cumtrapz(R**2/x,x)
    # phi[1:] = -1/x[1:]*cumulative_simpson(R**2,x=x) - cumulative_simpson(R**2/x,x=x)
    phi[0] = phi[1] #simps težave pri integraciji arraya dolžine nič, trapz pa ne!
    return phi

z=2 #št. elektronov
z0=z-5/16
def E_bohr(Z,n,E0=13.6058):
    '''energijski nivoji po Bohrovem modelu atoma,
    Z: atomsko (vrstno) št.'''
    return -E0 * Z**2/n**2
En=np.vectorize(E_bohr)

#  grid search
###########################################
# x=np.linspace(1e-9,10,200)
# R0 = lambda x: 2*z0**0.5*x*z0*np.exp(-z0*x) #zač. približek navodila
# phi = calc_potential(R0(x), x)

# y_zp = np.array([0,0.1]) #RK - lega,naklon
# # y_zpn = np.array([0,R0(x[0])])  #Numerov

# y_zp_inf=np.array([0, 0.1]) #neskončnost - integracija nazaj
# x1=x[::-1]

# epsilon_list= -np.linspace(0,10,500)
# R_list=[]
# Rn_list=[]
# R_nazaj_list=[]
# Rn_nazaj_list=[]
# for epsilon in epsilon_list:
#     # res = solve_ivp(DE, t_span=(x_min,x_max), y0=y_zp, t_eval=x, method='RK45',args=(epsilon,phi,z))
#     # R, dRdx = res.y #težava z phi!
#     R=numerov(x,y_zp,epsilon,phi)
#     Rn_list.append(R[-1])
#     R,_=runge_kutta_4(DE, x, y_zp, epsilon, phi,z)
#     R_list.append(R[-1])
#     R=numerov(x1,y_zp_inf,epsilon,phi)
#     Rn_nazaj_list.append(R[-1])
#     R,_=runge_kutta_4(DE, x1, y_zp_inf, epsilon, phi,z)
#     R_nazaj_list.append(R[-1])

# plt.plot(epsilon_list, np.abs(R_list),label='RK4 naprej')
# plt.plot(epsilon_list, np.abs(Rn_list),label='Numerov naprej')
# plt.plot(epsilon_list, np.abs(R_nazaj_list),label='RK4 nazaj')
# plt.plot(epsilon_list, np.abs(Rn_nazaj_list),label='Numerov nazaj')
# plt.legend()
# plt.yscale('log')
# plt.xlabel('$\\varepsilon$')
# plt.ylabel('$|R[-1]|$')
# plt.show()

# iteracije
#############################################
def last_f(epsilon, method, phi, N, x_max, y2=0, z=2, x_min=1e-9, R_return=False):
    y_zp_inf=np.array([1e-6, 0])
    # y_zp_infn=np.array([1e-6, 1e-6]) #numerov podani prvi 2 točki
    y_zp_0=np.array([0, 0.1])
    # y_zp_0n=np.array([0,y2])
    if method=='RK4 naprej':
        x=np.linspace(x_min,x_max, N)
        R,_ = runge_kutta_4(DE, x, y_zp_0, epsilon, phi, z)
    elif method=='RK4 nazaj':
        x=np.linspace(x_max,x_min, N)
        R,_ = runge_kutta_4(DE, x, y_zp_inf, epsilon, phi, z)
    elif method=='Numerov nazaj':
        x=np.linspace(x_max,x_min, N)
        R = numerov(x, y_zp_inf, epsilon, phi, z)
    elif method=='Numerov naprej':
        x=np.linspace(x_min,x_max, N)
        R = numerov(x, y_zp_0, epsilon, phi, z)
    if R_return:
        if x[1] < x[0]:
            R=R[::-1]
        return R
    else: return R[-1] #napaka na RP -> abs
def nicle(method, epsilon_0, epsilon_1, R0, phi, N=100, x_max=10, x_min=1e-9, z=2):
    y2=R0[1]
    try:
        sol = opt.root_scalar(last_f, x0=epsilon_0, x1=epsilon_1, method='secant', xtol=10**-10, args=(method,phi,N, x_max, y2, z), maxiter=200)
        epsilon = sol.root #bolj delikatni zač. približki
        print('E:', epsilon, 'secant', sol.iterations, sol.converged)
        if not sol.converged: raise ValueError('Sekantna metoda ne konvergira.')
    except ValueError:
        # sol = opt.root(last_f, x0=E0, method='broyden1', tol=10**-10, args=(method,phi,N,x_max,y2), options={'maxiter':1000})
        # sol = opt.root(last_f, x0=E0, method='hybr', tol=10**-10, args=(method,phi,N,x_max,y2), options={'maxiter':1000})
        sol = opt.root(last_f, x0=epsilon_0, method='lm', tol=10**-10, args=(method,phi,N,x_max,y2,z))
        epsilon = sol.x[0]
        print('E:', epsilon, 'nelin', sol.nfev, sol.success)
    return epsilon
def energija(R,x, epsilon, phi, E0=13.6058, z=2):
    '''skupna energija'''
    E = 2*E0 * simps(np.gradient(R,x)**2 - 2*R**2*z/x - phi*R**2, x)
    E -= epsilon*E0*(simps(R**2,x)-1) #normiranje
    return E
@njit
def count(y):
    '''štetje ničel'''
    n=1 #ena ničla na začetku
    for i in range(2,len(y)-1-1):
        if ((y[i]>=0 and y[i+1]<=0 and y[i]!=y[i+1]) or
            (y[i]<=0 and y[i+1]>=0 and y[i]!=y[i+1])):
            n+=1
    return n
def HF(method, N_iter, epsilon_0, R0, dx, x_max, z=2, x_min=1e-9, plot_iter=False, tol=1e-10):
    N=round((x_max-x_min)/dx)
    x=np.linspace(x_min,x_max,N)
    R=R0
    # R/=np.sqrt(simps(R**2,x)) #ZP ne normiramo
    phi = calc_potential(R0, x)
    en = energija(R0,x,epsilon_0, phi,z=z) #skupna E
    Rx=last_f(epsilon_0,method,phi,N,x_max,R0[1],R_return=True,z=z)
    Rx/=np.sqrt(simps(Rx[::-1]**2,x))

    # fig,ax=plt.subplots()
    # ax.plot(x, R, 'C0',label='$R_0(x)$')
    # ax.plot(x, Rx, 'k--', label='$R(x)$ pri $E_0$')
    # ax.set_xlabel('$x$')
    # ax.set_ylabel('$R_0(x)$', color='C0')
    # ax.tick_params(axis='y',labelcolor='C0')
    # ax.legend()
    # ax1=plt.twinx(ax)
    # ax1.plot(x, phi, 'C1', label='$\Phi_0(x)$')
    # ax1.set_ylabel('$\Phi_0(x)$', color='C1')
    # ax1.tick_params(axis='y',labelcolor='C1')
    # ax.set_title(f'zač. približek, $\\varepsilon_0={epsilon_0:<.3g}$, $E={en:.3f}$ eV')
    # alignYaxes([ax,ax1], [1,0])
    # plt.tight_layout()
    # plt.show()

    # epsilon_grid= np.linspace(0,10,50)
    # x1=np.linspace(x_min, x_max, 100) #za grid search

    epsilon_list = np.zeros(N_iter+1)
    epsilon_list[0]=epsilon_list[-1]=epsilon_0
    R_list=[R]
    energije=[en]
    epsilon_rel=[]
    energije_rel=[]
    epsilon=epsilon_0*1.001
    i=0
    e_relative=1
    while i<N_iter and tol<e_relative:
        phi=calc_potential(R,x)

        # grid search
        # R_grid=[]
        # for Ei in epsilon_grid:
        #     Ri=last_f(Ei,method,phi,N,x_max, y2=R[0])
        #     R_grid.append(Ri)
        # plt.plot(epsilon_grid, np.abs(R_grid),label=method)
        # plt.legend()
        # plt.yscale('log')
        # plt.xlabel('$\\varepsilon$')
        # plt.ylabel('$|R[-1]|$')
        # plt.title(f'iteracija {i+1}')
        # plt.show()

        epsilon = nicle(method, epsilon_list[i-1], epsilon, R, phi, N=N, z=z)
        epsilon_list[i+1]=epsilon

        R=last_f(epsilon,method,phi,N,x_max, y2=R[1], R_return=True, z=z)
        R /= np.sqrt(simps(R**2,x)) #np.sqrt(np.sum(R**2 *x)) #normalizacija
        R_list.append(R)
        en=energija(R,x,epsilon,phi,z=z)
        energije.append(en)
        e_relative = abs(energije[i]-energije[i-1])
        epsilon_rel.append(abs(epsilon_list[i]-epsilon_list[i-1]))
        energije_rel.append(e_relative)

        if plot_iter:
            fig,ax=plt.subplots()
            ax.plot(x, R, 'C0',label='$R_0(x)$')
            ax.set_xlabel('$x$')
            ax.set_ylabel('$R_0(x)$', color='C0')
            ax.tick_params(axis='y',labelcolor='C0')
            ax1=plt.twinx(ax)
            ax1.plot(x, np.abs(phi), 'C1', label='$\Phi_0(x)$')
            ax1.set_ylabel('$|\Phi_0(x)|$', color='C1')
            ax1.tick_params(axis='y',labelcolor='C1')
            ax1.set_yscale('log')
            ax.set_title(f'iteracija {i+1}: $\\varepsilon={epsilon:.3f}$, $E={en:.3f}$ eV')
            # alignYaxes([ax,ax1])
            plt.show()
        i+=1
    if i>=N_iter: print(f'Max iteracije {N_iter} dosežene!')
    else:
        epsilon_list=epsilon_list[:i+1] #trim zeros
        print(f'Toleranca {tol:.0E} dosežena ({e_relative}) v {i} iteracijah.')
    nz = count(R)
    return epsilon_list, R_list, energije, epsilon_rel, energije_rel, nz, i

# začetni približki
#############################################################################
method='Numerov nazaj'
N_iter = 20
x_min, x_max, dx = 1e-9, 10, 0.01
epsilon_0 = -4 #strelski param., lagrangeev multiplikator
R0 = lambda x: 2*z0**0.5*x*z0*np.exp(-z0*x) #zač. približek

x=np.linspace(x_min,x_max,round((x_max-x_min)/dx))
R0 = R0(x)
epsilon_list, R_list, energije, epsilon_rel, energije_rel, nz,i = HF(method, N_iter, epsilon_0, R0, dx, x_max)

plt.plot(epsilon_list, '.-',label='$\\varepsilon$')
plt.plot(energije, '.-',label='$E$')
plt.legend()
plt.xlabel('iteracije')
plt.ylabel('$E$ [eV], $\\varepsilon$ [/]')
plt.title(method)
plt.show()

plt.plot(np.arange(1,len(epsilon_list),1), epsilon_rel,'.-',label='$\\varepsilon$')
plt.plot(np.arange(1,len(epsilon_list),1), energije_rel,'.-',label='$E$')
plt.yscale('log')
plt.xlabel('iteracije')
plt.ylabel('$|E_i-E_{i-1}|$, $|\\varepsilon_i-\\varepsilon_{i-1}|$ [eV]')
plt.legend()
plt.show()

cmap=plt.get_cmap('viridis')
norm=Normalize(1,len(epsilon_list)-1)
fig, ax = plt.subplots()
for i in range(1,len(R_list)):
    # plt.plot(x, R_list[i], c=cmap(norm(i)),alpha=0.6)
    plt.plot(x, np.abs(R_list[i]-R_list[i-1]), c=cmap(norm(i)),alpha=0.9)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
plt.colorbar(sm,ax=ax, label='iteracije')
plt.xlabel('$x$')
# plt.ylabel('$R(x)$')
plt.ylabel('$|R(x)-R_0(x)|$')
# plt.savefig(f'./slike/konvergenca_R.pdf')
plt.show()

# # nov zač. približek
# R0 = np.full_like(x, 1./x_max) #konstanta - normirano
# epsilon_list1, R_list1, energije1, epsilon_rel1, energije_rel1, nz1,i1 = HF(method, N_iter, epsilon_0, R0, dx, x_max)

# R0 = np.random.random(len(x))
# epsilon_list2, R_list2, energije2, epsilon_rel2, energije_rel2, nz2,i2 = HF(method, N_iter, epsilon_0, R0, dx, x_max)

# print('končne vrednosti epsilon:', epsilon_list[-1], epsilon_list1[-1], epsilon_list2[-1])
# print('končne vrednosti energije:', energije[-1],energije1[-1],energije2[-1])
# print('končne vrednosti iteracije:', i,i1,i2)
# ax = plt.figure().gca()
# plt.plot(np.arange(1,len(epsilon_list),1), epsilon_rel,'C0.-',label='$R_0(x)=2 \sqrt{Z^*} x Z^* e^{-Z^* x}$')
# plt.plot(np.arange(1,len(epsilon_list),1), energije_rel,'C0.--')

# plt.plot(np.arange(1,len(epsilon_list1),1), epsilon_rel1,'C1.-',label='$R_0(x)=konst.={1./x_max}$')
# plt.plot(np.arange(1,len(epsilon_list1),1), energije_rel1,'C1.--')

# plt.plot(np.arange(1,len(epsilon_list2),1), epsilon_rel2,'C2.-',label='$R_0(x)=$Random')
# plt.plot(np.arange(1,len(epsilon_list2),1), energije_rel2,'C2.--')

# ax.xaxis.set_major_locator(MaxNLocator(integer=True))
# plt.title(f'polna črta: skupna $E$, prekinjena: $\\varepsilon$')
# plt.yscale('log')
# plt.xlabel('iteracije')
# plt.ylabel('$|E_i-E_{i-1}|$, $|\\varepsilon_i-\\varepsilon_{i-1}|$ [eV]')
# plt.legend()
# # plt.savefig(f'./slike/konvergenca_E_epsilon_ZP.pdf')
# plt.show()

# plt.plot(x, R_list[-1])
# plt.plot(x, R_list1[-1])
# plt.plot(x, R_list2[-1])
# plt.show()

# različni integratorji
#######################################################
# N_iter = 50
# x_min, x_max, dx = 1e-9, 10, 0.01
# epsilon_0 = -5 #strelski param., lagrangeev multiplikator
# R0 = lambda x: 2*z0**0.5*x*z0*np.exp(-z0*x) #zač. približek
# x=np.linspace(x_min,x_max,round((x_max-x_min)/dx))
# R0 = R0(x)

# method='Numerov nazaj'
# t1=timer()
# epsilon_list, R_list, energije, epsilon_rel, energije_rel, nz,i = HF(method, N_iter, epsilon_0, R0, dx, x_max)
# t2=timer()
# time0=t2-t1

# method='Numerov naprej'
# t1=timer()
# epsilon_list1, R_list1, energije1, epsilon_rel1, energije_rel1, nz1,i1 = HF(method, N_iter, epsilon_0, R0, dx, x_max)
# t2=timer()
# time1=t2-t1

# method='RK4 nazaj'
# t1=timer()
# epsilon_list2, R_list2, energije2, epsilon_rel2, energije_rel2, nz2,i2 = HF(method, N_iter, epsilon_0, R0, dx, x_max)
# t2=timer()
# time2=t2-t1

# method='RK4 naprej'
# t1=timer()
# epsilon_list3, R_list3, energije3, epsilon_rel3, energije_rel3, nz3,i3 = HF(method, N_iter, epsilon_0, R0, dx, x_max)
# t2=timer()
# time3=t2-t1

# ax = plt.figure().gca()
# plt.plot(epsilon_rel,'C0.-',label='Numerov nazaj')
# plt.plot(energije_rel,'C0.--')
# plt.plot(epsilon_rel1,'C1.-',label='Numerov naprej')
# plt.plot(energije_rel1,'C1.--')
# plt.plot(epsilon_rel2,'C2.-',label='RK4 nazaj')
# plt.plot(energije_rel2,'C2.--')
# plt.plot(epsilon_rel3,'C3.-',label='RK4 naprej')
# plt.plot(energije_rel3,'C3.--')

# ax.xaxis.set_major_locator(MaxNLocator(integer=True))
# plt.title(f'polna črta: skupna $E$, prekinjena: $\\varepsilon$')
# plt.yscale('log')
# plt.xlabel('iteracije')
# plt.ylabel('$|E_i-E_{i-1}|$, $|\\varepsilon_i-\\varepsilon_{i-1}|$ [eV]')
# plt.legend()
# # plt.savefig(f'./slike/konv_E_epsilon_integratorji.pdf')
# plt.show()

# epsilon_data=[epsilon_list[-1],epsilon_list1[-1],epsilon_list2[-1],epsilon_list3[-1]]
# energije_data=[energije[-1],energije1[-1],energije2[-1],energije3[-1]]
# iters=[i,i1,i2,i3]
# n_nicle=[nz,nz1,nz2,nz3]
# time=[time0,time1,time2,time3]
# R_data=[R_list[-1],R_list1[-1],R_list2[-1],R_list3[-1]]
# metode=['Numerov nazaj', 'Numerov naprej', 'RK4 nazaj', 'RK4 naprej']

# ax = plt.figure().gca()
# ax.plot(iters, '+', label='iteracije', markersize=10)
# ax.plot(n_nicle, '+', label='ničle', markersize=10)
# ax.yaxis.set_major_locator(MaxNLocator(integer=True))
# ax.set_xticks([0,1,2,3], labels=metode)
# ax.set_ylabel('št. iteracij / ničel')
# ax1=plt.twinx(ax)
# ax1.plot( time, 'g+', markersize=10)
# ax1.tick_params(axis='y',labelcolor='g')
# ax1.set_ylabel('$t$ [s]', color='g')
# ax.legend(loc='upper left')
# alignYaxes([ax,ax1],[0,0])
# # plt.savefig(f'./slike/integratorji_iter.pdf')
# plt.show()

# for i in range(len(R_data)):
#     plt.plot(x, R_data[i], label=f'{metode[i]}: $\\varepsilon={epsilon_data[i]:.3g}$, $E={energije_data[i]:.3g}$')
# plt.legend()
# plt.xlabel('$x$')
# plt.ylabel('$R(x)$')
# # plt.savefig(f'./slike/integratorji_R.pdf')
# plt.show()

# odvisnost od dx
#########################################################
Z=2
method='Numerov nazaj'
N_iter = 100
x_min, x_max, dx = 1e-9, 10, 0.01
epsilon_0=-5
R0 = lambda x: 2*z0**0.5*x*z0*np.exp(-z0*x) #zač. približek

epsilon_vals=[]
e_vals=[]
n_nicle=[]
iters=[]
R_data=[]
R0_data=[]
time=[]
e_end_vals=[]
epsilon_end=[]
dx_data=np.linspace(0.001,0.5,20)
for dx in dx_data:
    x=np.linspace(x_min,x_max,round((x_max-x_min)/dx))
    R0i = R0(x)
    t1=timer()
    epsilon_list, R_list, energije, epsilon_rel, energije_rel, nz,i = HF(method, N_iter, epsilon_0, R0i, dx, x_max, Z)
    t2=timer()
    epsilon_vals.append(epsilon_rel)
    e_vals.append(energije_rel)
    R_data.append(R_list[-1])
    n_nicle.append(nz)
    iters.append(i)
    time.append(t2-t1)
    e_end_vals.append(energije[-1])
    epsilon_end.append(epsilon_list[-1])

plt.plot(dx_data, np.abs(e_end_vals),'.-',label='$E$')
plt.plot(dx_data, np.abs(epsilon_end),'.-',label='$\\varepsilon$')
plt.xscale('linear')
plt.yscale('log')
plt.legend()
plt.xlabel('$\Delta x$')
plt.ylabel('$E$, $\\varepsilon$ [eV]')
# plt.savefig('./slike/konv_dx_E.pdf')
plt.show()

cmap=plt.get_cmap('viridis')
norm=Normalize(min(dx_data),max(dx_data))
fig, ax = plt.subplots()
for i in range(len(dx_data)):
    plt.plot(e_vals[i], c=cmap(norm(dx_data[i])),alpha=0.9)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
plt.colorbar(sm,ax=ax, label='$\Delta x$')
plt.xlabel('$x$')
# plt.ylabel('$R(x)$')
plt.ylabel('$|E_i-E_{i-1}|$')
plt.yscale('log')
# plt.savefig(f'./slike/konv_dx_E_rel.pdf')
plt.show()

# več ničel -> višja vzbujena stanja
########################################################
# Z=2
# method='Numerov nazaj'
# method1='Numerov naprej'
# N_iter = 50
# x_min, x_max, dx = 1e-9, 50, 0.05 #0.01 bolje, a zelo počasi; 0.1 grobo?
# R0 = lambda x: 2*z0**0.5*x*z0*np.exp(-z0*x) #zač. približek
# x=np.linspace(x_min,x_max,round((x_max-x_min)/dx))
# R0 = R0(x)
# epsilon_vals=[]
# e_vals=[]
# n_nicle=[]
# iters=[]
# R_data=[]
# time=[]
# # epsilon_vals1=[]
# # e_vals1=[]
# # n_nicle1=[]
# # iters1=[]
# # R_data1=[]
# # time1=[]
# # epsilon_data=-np.logspace(-13,17,20, base=1.1)
# epsilon_data=-np.linspace(1e-5,5,20)
# # input(epsilon_data)
# for epsilon_0 in epsilon_data:
#     t1=timer()
#     epsilon_list, R_list, energije, epsilon_rel, energije_rel, nz,i = HF(method, N_iter, epsilon_0, R0, dx, x_max, Z)
#     t2=timer()
#     epsilon_vals.append(epsilon_list[-1])
#     e_vals.append(energije[-1])
#     if nz not in n_nicle:
#         R_data.append(R_list[-1])
#         # plt.plot(x,R_list[-1])
#         # plt.title(f'$\\varepsilon={epsilon_list[-1]:.4g}$, $nz={nz}$')
#         # plt.show()
#     n_nicle.append(nz)
#     iters.append(i)
#     time.append(t2-t1)

#     # #metoda1
#     # t1=timer()
#     # epsilon_list, R_list, energije, epsilon_rel, energije_rel, nz,i = HF(method1, N_iter, epsilon_0, R0, dx, x_max, Z)
#     # t2=timer()
#     # epsilon_vals1.append(epsilon_list[-1])
#     # e_vals1.append(energije[-1])
#     # if nz not in n_nicle1:
#     #     R_data1.append(R_list[-1])
#     #     # plt.plot(x,R_list[-1])
#     #     # plt.title(f'$\\varepsilon={epsilon_list[-1]:.4g}$, $nz={nz}$')
#     #     # plt.show()
#     # n_nicle1.append(nz)
#     # iters1.append(i)
#     # time1.append(t2-t1)

# e_vals=np.array(e_vals,dtype=np.float32) #izogib numerične napake
# e_nivoji=sorted(list(set(e_vals))) #diskretni energijski nivoji
# epsilon_vals=np.array(epsilon_vals,dtype=np.float32) #izogib numerične napake
# epsilon_vals=sorted(list(set(epsilon_vals)))
# nicle_sort=sorted(list(set(n_nicle)))

# # e_vals1=np.array(e_vals1,dtype=np.float32) #izogib numerične napake
# # e_nivoji1=sorted(list(set(e_vals1))) #diskretni energijski nivoji
# # epsilon_vals1=np.array(epsilon_vals1,dtype=np.float32) #izogib numerične napake
# # epsilon_vals1=sorted(list(set(epsilon_vals1)))
# # nicle_sort1=sorted(list(set(n_nicle1)))

# for i in range(len(e_nivoji)):
#     plt.axhline(e_nivoji[i], c='k', ls='--')
#     plt.text(0.3,e_nivoji[i],f'$n={nicle_sort[i]-1}$')
#     plt.axvline(epsilon_vals[i], c='b', ls='--')
# plt.plot(epsilon_data, e_vals, 'rx')
# plt.xlabel('$\\varepsilon_0$ [eV]')
# plt.ylabel('skupna $E$ [eV]')
# plt.title(method)
# plt.tight_layout()
# # plt.savefig(f'./slike/vecji_n_energije_Z{Z}_{method}.pdf')
# plt.show()
# # for i in range(len(e_nivoji1)):
# #     plt.axhline(e_nivoji1[i], c='k', ls='--')
# #     plt.text(0.3,e_nivoji1[i],f'$n={nicle_sort1[i]-1}$')
# #     plt.axvline(epsilon_vals1[i], c='b', ls='--')
# # plt.plot(epsilon_data, e_vals1, 'rx')
# # plt.xlabel('$\\varepsilon_0$ [eV]')
# # plt.ylabel('skupna $E$ [eV]')
# # plt.title(method1)
# # plt.tight_layout()
# # # plt.savefig(f'./slike/vecji_n_energije_Z{Z}_{method1}.pdf')
# # plt.show()

# ax = plt.figure().gca()
# ax.plot(epsilon_data, iters, 'C0.-', label='iteracije')
# ax.plot(epsilon_data, n_nicle, 'C1.-', label='ničle')
# # ax.plot(epsilon_data, iters1, 'C0*--')
# # ax.plot(epsilon_data, n_nicle1, 'C1*--')
# ax.yaxis.set_major_locator(MaxNLocator(integer=True))
# ax.set_xlabel('$\\varepsilon_0$ [eV]')
# ax.set_ylabel('št. iteracij / ničel')
# ax1=plt.twinx(ax)
# ax1.plot(epsilon_data, time, 'g.-')
# ax1.tick_params(axis='y',labelcolor='g')
# ax1.set_ylabel('$t$ [s]', color='g')
# ax.legend(loc='upper left')
# # plt.title(f'polna črta: Numerov naprej, prekinjena: nazaj')
# alignYaxes([ax,ax1],[0,0])
# # plt.savefig(f'./slike/vecji_n_iters_Z{Z}.pdf')
# plt.show()

# c=['C0','C1','C2','C3','C4','C5','C6','C7']
# for i in range(len(R_data)):
#     plt.plot(x, R_data[i], label=f'{nicle_sort[-i-1]}',c=c[i])
# plt.legend(title='lastne funkcije $n$')
# plt.xlabel('$x$')
# plt.ylabel('$R(x)$')
# # plt.savefig(f'./slike/vecji_n_R_{method}_Z{Z}.pdf')
# plt.show()
# # for i in range(len(R_data1)):
# #     plt.plot(x, R_data1[i], label=f'{nicle_sort1[-i-1]}',c=c[i])
# # plt.legend(title='lastne funkcije $n$')
# # plt.xlabel('$x$')
# # plt.ylabel('$R(x)$')
# # # plt.savefig(f'./slike/vecji_n_R_{method1}.pdf')
# # plt.show()

# ax = plt.figure().gca()
# # n=np.linspace(1,len(e_nivoji), len(e_nivoji))
# plt.plot(nicle_sort, e_nivoji, '.-', label='Hartree-Fock: Numerov nazaj')
# # plt.plot(nicle_sort1, e_nivoji1, '.-', label='Hartree-Fock: Numerov naprej')
# plt.plot(nicle_sort, E_bohr(Z=Z, n=np.array(nicle_sort)), '.-', label='Bohr')
# ax.xaxis.set_major_locator(MaxNLocator(integer=True))
# plt.xlabel('$n$')
# plt.ylabel('$E(n)$ [eV]')
# plt.legend()
# # plt.savefig(f'./slike/primerjava_Bohr_Z{Z}.pdf')
# plt.show()

# različno št e- -> Z
#########################################################
method='Numerov nazaj'
N_iter = 100
x_min, x_max, dx = 1e-9, 10, 0.01
R0 = lambda x: 2*z0**0.5*x*z0*np.exp(-z0*x) #zač. približek
x=np.linspace(x_min,x_max,round((x_max-x_min)/dx))
R0 = R0(x)

# Z=20
# epsilon_0=-Z**2
# epsilon_list, R_list, energije, epsilon_rel, energije_rel, nz,i = HF(method, N_iter, epsilon_0, R0, dx, x_max, Z)

# plt.plot(x,R_list[-1])
# plt.show()
# plt.plot(np.abs(epsilon_list),'.-')
# plt.plot(np.abs(energije),'.-')
# plt.yscale('log')
# plt.show()
# plt.plot(epsilon_rel,'.-')
# plt.plot(energije_rel,'.-')
# plt.yscale('log')
# plt.show()

############ zanka
epsilon_vals=[]
e_vals=[]
n_nicle=[]
iters=[]
R_data=[]
time=[]
Z_list=np.linspace(1,10,20)
# Z_list=[1,2,3,4,5,10] #za prikaz
for Z in Z_list:
    epsilon_0=-Z**2
    t1=timer()
    epsilon_list, R_list, energije, epsilon_rel, energije_rel, nz,i = HF(method, N_iter, epsilon_0, R0, dx, x_max, Z)
    t2=timer()
    epsilon_vals.append(epsilon_list[-1])
    e_vals.append(energije[-1])
    R_data.append(R_list[-1])
    # if nz not in n_nicle:
        # plt.plot(x,R_list[-1])
        # plt.title(f'$Z={Z}$, $\\varepsilon={epsilon_list[-1]:.4g}$, $nz={nz}$')
        # plt.show()
    n_nicle.append(nz)
    iters.append(i)
    time.append(t2-t1)
print(set(n_nicle))
# print(method, e_vals[0], e_vals[1], e_vals[2])

e_vals=np.array(e_vals,dtype=np.float32) #izogib numerične napake
e_nivoji=sorted(list(set(e_vals))) #diskretni energijski nivoji
epsilon_vals=np.array(epsilon_vals,dtype=np.float32) #izogib numerične napake
epsilon_vals=sorted(list(set(epsilon_vals)))
nicle_sort=sorted(list(set(n_nicle)))

plt.plot(Z_list, np.abs(e_vals), 'k.-', label='izračun')
plt.plot([1,2,3],np.abs([-14.34,-78.88,-198.04]), 'rx', label='eksperiment')
plt.legend()
plt.yscale('log')
plt.xlabel('$Z$')
plt.ylabel('skupna $|E|$ [eV]')
plt.tight_layout()
# plt.savefig(f'./slike/Z_energija.pdf')
plt.show()

ax = plt.figure().gca()
ax.plot(Z_list, iters, '.-', label='itearcije')
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
ax.set_xlabel('$Z$')
ax.set_ylabel('št. iteracij')
ax1=plt.twinx(ax)
ax1.plot(Z_list, time, 'g.-')
ax1.tick_params(axis='y',labelcolor='g')
ax1.set_ylabel('$t$ [s]', color='g')
alignYaxes([ax,ax1],[0,0])
# plt.savefig(f'./slike/Z_iteracije.pdf')
plt.show()

for i in range(len(R_data)):
    plt.plot(x, R_data[i], label=f'{Z_list[i]}')
plt.legend(title='$Z$')
plt.xlabel('$x$')
plt.ylabel('$R(x)$')
plt.title(method)
# plt.savefig(f'./slike/R_funk_Z_naprej.pdf')
plt.show()