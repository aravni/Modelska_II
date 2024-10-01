# -*- coding: utf-8 -*-
# 3. naloga lastne energije schroedingerjeve en.
# 1. del: radialne del en. za vodikov atom s Coulombskim potencialom
# ---- povzeto po sprotnem reševanju iz datoteke 3.nal.ipynb ----

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.integrate import solve_ivp
from scipy import optimize as opt
from scipy.integrate import simps
from scipy.special import genlaguerre
from timeit import default_timer as timer
from numba import njit, jit
from matplotlib.colors import ListedColormap, BoundaryNorm

plt.style.use('aljaz')
@njit
def DE(x,y, e,l=0):
    '''integracija en. za vodikov atom'''
    R,z = y.reshape((2,len(y)//2))
    dRdx = z
    dzdx = (-2/x + (l*(l+1))/x**2 - e) * R
    return np.append(dRdx, dzdx)
@njit
def k_fun(x,e,l):
    '''funkcija k²(x) za metodo Numerov, Schroedingerjeva en. vodikov atom'''
    return 2/x - l*(l+1)/x**2 + e
@njit
def R_apr(x,koef):
    '''polinomska aproksimacija rešitve Schroedingerjeve en.'''
    y=0
    for i in range(len(koef)):
        y += koef[i]*x**(i+1)
    return y
@njit
def numerov(x, y0, e, l=0):
    '''reševanje DE oblike [d²/dx² + k²(x)] y(x) = 0,
    funkcija k(x) naj bo že kvadrirana. Pazi na časovni korak!'''
    h = x[1]-x[0]
    y = np.zeros_like(x)

    y[0] = y0[0] #ZP
    y[1] = y0[0] + h*y0[1] #lin. razvoj -> OK

    # N=10 #red polinoma za aproksimacijo rešitve -> ne dela (slabo) za podlupine
    # a = np.zeros(N) #koef. polinomske aproksimacije rešitve
    # if l==0:
    #     a[0] = 1
    #     a[1] = -a[0]
    #     for i in range(2,N):
    #         a[i] = -1/(i*(i-1)) * (2*a[i-1]+e*a[i-2])
    # y[0] = R_apr(x[0], a)
    # y[1] = R_apr(x[1], a)

    for i in range(1,len(x)-1):
        a_2 = 1 - h**2*5/12*k_fun(x[i],e,l)
        a_1 = 1 + h**2*1/12*k_fun(x[i-1],e,l)
        a_3 = 1 + h**2*1/12*k_fun(x[i+1],e,l)
        y[i+1] = (2*a_2*y[i] - a_1*y[i-1]) / a_3
    return y
@njit
def runge_kutta_4(f, t, y0, e, l):
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
        k0 = h*f(ti, y[i], e, l)
        k1 = h*f(ti + h/2.0, y[i] + k0/2.0, e, l)
        k2 = h*f(ti + h/2.0, y[i] + k1/2.0, e, l)
        k3 = h*f(ti + h, y[i] + k2, e, l)
        y[i+1] = y[i] + (k0 + 2.0*k1 + 2.0*k2 + k3)/6.0
        i+=1
    return np.transpose(y)

# točne (nenormalizirane) rešitve - last. funk.
# prva št. = n (glavno kvantno število, lupina), n > 0
# drugo št = l (tirno kvantno število-vrtilna kol., podlupina), 0 <= l <= n-1
def R10(x):
    y=x*np.e**(-x)
    dx = abs(x[1] - x[0])
    norm = np.sqrt(np.sum(np.abs(y)**2 * dx)) 
    y /= norm
    return y
def R20(x):
    y=x*(1-x/2)*np.e**(-x/2)
    dx = abs(x[1] - x[0])
    norm = np.sqrt(np.sum(np.abs(y)**2 * dx)) 
    y /= norm
    return y
def R21(x):
    y=x**2*np.e**(-x/2)
    dx = abs(x[1] - x[0])
    norm = np.sqrt(np.sum(np.abs(y)**2 * dx)) 
    y /= norm
    return y
def R30(x):
    y=x*np.exp(-x/3)*(1-2/3*x+2/27*x**2)
    dx = abs(x[1] - x[0])
    norm = np.sqrt(np.sum(np.abs(y)**2 * dx)) 
    y /= norm
    return y
def R31(x):
    y=x*x*np.exp(-x/3)*(1-1/6*x)
    dx = abs(x[1] - x[0])
    norm = np.sqrt(np.sum(np.abs(y)**2 * dx)) 
    y /= norm
    return y
def R32(x):
    y=x**3*np.exp(-x/3)
    dx = abs(x[1] - x[0])
    norm = np.sqrt(np.sum(np.abs(y)**2 * dx)) 
    y /= norm
    return y
def Rnl(r, n, l, a0=1): #NI OK
    '''Radial wavefunction of hydrogen atom for quantum numbers n and l'''
    rho = 2 * r / (n * a0)
    norm = np.sqrt((2 / (n * a0))**3 * np.math.factorial(n - l - 1) / (2 * n * np.math.factorial(n + l)))
    laguerre = genlaguerre(n - l - 1, 2 * l + 1)(rho)
    return norm * np.exp(-rho / 2) * (rho**l) * laguerre

# poskusno streljanje
####################################################################
n=3 #lupina
l=2 #podlupina
x_max = 10*n**1.5 #x_max prilagodi željeni lupini!!
y_0_inf = 1e-4 #ZP v neskončosti
x = np.linspace(1e-9,x_max,500)# radialna koordinata
dx = x[1] - x[0]  # Spacing between x points
x_min=x[0] #singularno v x=0

y_zp = np.array([0, 1]) #začetna lega, naklon -> naklon NE vpliva na obliko rešitve -> normalizacija
e0 = -1/n**2 #analitična lastna energija -> shooting parameter

# integracija naprej: občutlivo na prevelik x_max, ker integrator pobezlja v neskončnosti
res = solve_ivp(DE, t_span=(x_min,x_max), y0=y_zp, t_eval=x, method='RK45',args=(e0,l))
R1, dR1dx = res.y
x1 = res.t
norm = np.sqrt(np.sum(np.abs(R1)**2 * dx))  # Approximate the integral
R1 /= norm  # Normalize the wavefunction

# integracija nazaj: občutljivo na premajhen x_max, da funk ni dovolj blizu nič
y_zp = np.array([y_0_inf, 0]) #zp v neskončnosti!
res = solve_ivp(DE, t_span=(x_max,x_min), y0=y_zp, t_eval=x[::-1], method='RK45',args=(e0,l))
R2, dR2dx = res.y
x2 = res.t
if R2[-2]-R2[-1] < 0: # naklon na zač. želimo 1
    R2/=-1
norm = np.sqrt(np.sum(np.abs(R2)**2 * dx))  # Approximate the integral -> OK
R2 /= norm  # Normalize the wavefunction
# print(np.trapz(R2**2,x2)+1,simps(R2**2,x2)+1)

# obojestransko
if (n-l)%2==0: y_0_inf=-y_0_inf #sodo št ničel, konvergenca od spodaj
y_zp = np.array([0, 1])
res = solve_ivp(DE, t_span=(x_min,x[len(x)//2]), y0=y_zp, method='RK45', args=(e0,l), t_eval=x[:len(x)//2+1]) #naprej
y1, dy1dx = res.y
x11=res.t 
y_zp = np.array([y_0_inf, 0]) #zp v neskončnosti
res = solve_ivp(DE, t_span=(x_max,x[len(x)//2]), y0=y_zp, method='RK45', args=(e0,l), t_eval=x[:len(x)//2-1:-1]) #nazaj
# res = solve_ivp(DE, t_span=(x_max,x_min), y0=y_zp, t_eval=x[::-1], method='RK45',args=(e0,l))
y2, dy2dx = res.y
x22=res.t
scale_factor = dy1dx[-1] / dy2dx[-1] #scale da je naklon enak
y2 = y2 * scale_factor
dy2dx *= scale_factor 
R3 = np.append(y1,y2[::-1]) #združitev
dR3dx = np.append(dy1dx,dy2dx[::-1]) #združitev
x3=np.append(x11,x22[::-1])
norm = np.sqrt(np.sum(np.abs(R3)**2 * dx))  # Approximate the integral
R3 /= norm  # Normalize the wavefunction

# NUMEROV
y_zp = np.array([0, 1])
# e0=-0.827895874989988 #n=1?
R4 = numerov(x,y_zp,e0,l)
x4=x
norm = np.sqrt(np.sum(np.abs(R4)**2 * dx))  # Approximate the integral
R4 /= norm  # Normalize the wavefunction

# # plt.plot(x, R10(x), label='R10')
# # plt.plot(x, R20(x), label='R20')
# # plt.plot(x, R21(x), label='R21')
# # plt.plot(x, R32(x), label='R32')
# plt.plot(x1, R1, label=f'naprej')
# plt.plot(x2, R2, label=f'nazaj')
# plt.plot(x3, R3, label=f'obojestransko')
# plt.plot(x4, R4, label=f'Numerov')
# plt.title(f'$e=-1/{n}^2$, $l={l}$')
# plt.legend()
# plt.show()

# plt.plot(x1, dR1dx, label=f'naprej')  #odvod = naklon
# plt.plot(x2, dR2dx, label=f'nazaj')
# plt.plot(x3, dR3dx, label=f'obojstransko')
# plt.title(f'$e=-1/{n}^2$, $l={l}$')
# plt.legend()
# plt.show()

# grid search
#####################################################################
# def napaka_energija(e, x_max, method='nazaj', l=0, y_0_inf=1e-4):
#     '''vrne poz. rezultat za prikaz v log skali'''
#     y_zp = np.array([0, 1]) #zp pri x=0
#     if (n-l)%2==0: y_0_inf=-y_0_inf #sodo št ničel, konvergenca od spodaj
#     if method =='nazaj':
#         y_zp = np.array([y_0_inf, 0]) #zp v neskončnosti
#         y, dydx = solve_ivp(DE, t_span=(x_max,x_min), y0=y_zp, method='RK45', args=(e,l)).y #nazaj
#     elif method == 'naprej':
#         y, dydx = solve_ivp(DE, t_span=(x_min,x_max), y0=y_zp, method='RK45', args=(e,l)).y #naprej
#     elif method == 'numerov naprej':
#         x = np.linspace(x_min,x_max, 500)
#         y = numerov(x, y_zp, e)
#     elif method == 'numerov nazaj':
#         x = np.linspace(x_max,x_min, 500)
#         y_zp = np.array([y_0_inf, 0]) #zp v neskončnosti
#         y = numerov(x, y_zp, e)
#     elif method == 'obojestransko':
#         y1, dy1dx = solve_ivp(DE, t_span=(x_min,x_max/2), y0=y_zp, method='RK45', args=(e,l)).y #naprej
#         y_zp = np.array([y_0_inf, 0]) #zp v neskončnosti
#         y2, dy2dx = solve_ivp(DE, t_span=(x_max,x_max/2), y0=y_zp, method='RK45', args=(e,l)).y #nazaj
#         scale_factor = dy1dx[-1] / dy2dx[-1] #scale da je naklon enak
#         y2 = y2 * scale_factor
#         y = [y1[-1] - y2[-1]] #napaka zdužitve
#     return abs(y[-1]) #ok,  abs(y[-1]+ dydx[-1]) ni ok
#     # return (y[-1])**2 #ok,  (y[-1]+ dydx[-1])**2 ni ok
# # e = np.linspace(-1.1e0, -1e-6, 1000)
# e = -np.logspace(0.05, -2.5, 1000)

# t1=timer()
# R=np.zeros_like(e)
# for i,e_i in enumerate(e):
#     R[i] = napaka_energija(e_i, x_max=10*(-1/e_i)**(0.5), method='naprej')
# t2=timer()
# # napaka_energija_vec = np.vectorize(napaka_energija) #ni hitreje?
# # R = napaka_energija_vec(e,x_max=(-1/e_i)**0.5, method='naprej')
# # t1=timer()
# print(f'naprej t_sim = {abs(t2-t1):.1f}')
# R1=np.zeros_like(e)
# for i,e_i in enumerate(e):
#     R1[i] = napaka_energija(e_i, x_max=10*(-1/e_i)**(0.5), method='nazaj')
# t3=timer()
# print(f'nazaj t_sim = {(t3-t2):.1f}')
# R2=np.zeros_like(e)
# for i,e_i in enumerate(e):
#     R2[i] = napaka_energija(e_i, x_max=10*(-1/e_i)**(0.5), method='obojestransko')
# t4=timer()
# print(f'obojestransko t_sim = {(t4-t3):.1f}')
# Rn=np.zeros_like(e)
# for i,e_i in enumerate(e):
#     Rn[i] = napaka_energija(e_i, x_max=10*(-1/e_i)**(0.5), method='numerov naprej')
# print(f'numerov naprej t_sim = {(timer()-t4):.1f}')
# Rnn=np.zeros_like(e)
# for i,e_i in enumerate(e):
#     Rnn[i] = napaka_energija(e_i, x_max=10*(-1/e_i)**(0.5), method='numerov nazaj')
# print(f'numerov nazaj t_sim = {(timer()-t4):.1f}')

# plt.plot(-e,R,'-',label='RK naprej')
# plt.plot(-e,R1,'-',label='RK nazaj')
# plt.plot(-e,R2,'-',label='obojestransko')
# plt.plot(-e,Rn,'--',label='Numerov naprej')
# plt.plot(-e,Rnn,'--',label='Numerov nazaj')
# for ni in range(1,int((1/min(-e))**0.5)):
#     plt.axvline(1/ni**2, c='k', ls=':')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('$-E$')
# plt.ylabel('$|y_{i=-1}|$')
# plt.legend()
# # plt.savefig(f'./slike/grid_search.pdf')
# plt.show()

# iskanje za ZP
#############################################################
def napaka_energija(e, x_max, l=0, method='nazaj', N=100,  y_0_inf=1e-4, x_min=1e-9):
    '''za iskanje ničle'''
    y_zp = np.array([0, 1]) #zp pri x=0
    if (n-l)%2==0: y_0_inf=-y_0_inf #sodo št ničel, konvergenca od spodaj
    if method =='nazaj':
        y_zp = np.array([y_0_inf, 0]) #zp v neskončnosti
        y, _ = solve_ivp(DE, t_span=(x_max,x_min), y0=y_zp, method='RK45', args=(e,l)).y #nazaj
    elif method =='naprej-RK':
        x=np.linspace(x_min,x_max,N)
        y, _ = runge_kutta_4(DE, x, y_zp, e, l)
    elif method == 'naprej':
        y, _ = solve_ivp(DE, t_span=(x_min,x_max), y0=y_zp, method='RK45', args=(e,l)).y #naprej
    elif method == 'numerov':
        x=np.linspace(x_min,x_max,N)
        y = numerov(x, y_zp, e, l)
    elif method == 'numerov nazaj':
        x=np.linspace(x_max,x_min,N)
        y = numerov(x, y_zp, e, l)
    elif method == 'obojestransko':
        y1, dy1 = solve_ivp(DE, t_span=(x_min,x_max/2), y0=y_zp, method='RK45', args=(e,l)).y #naprej
        y_zp = np.array([y_0_inf, 0]) #zp v neskončnosti
        y2, dy2 = solve_ivp(DE, t_span=(x_max,x_max/2), y0=y_zp, method='RK45', args=(e,l)).y #nazaj
        scale_factor = dy1[-1] / dy2[-1] #scale da je naklon enak
        y2 = y2 * scale_factor
        y = [y1[-1]-y2[-1]] #napaka združitve
    else: print('napaka: ', method)
    return y[-1]
@njit
def count(y):
    n=1 #ena ničla na začetku, ena ničla odvoda na koncu
    m=max([np.max(y),-np.min(y)])
    for i in range(2,len(y)-1-1):
        if ((y[i]>=0 and y[i+1]<=0 and not np.isclose(y[i],y[i+1],atol=m*1e-9)) or
            (y[i]<=0 and y[i+1]>=0 and not np.isclose(y[i],y[i+1],atol=m*1e-9))):
            n+=1
    return n
def last_f(method, x_min,x_max, e0, l, N=100, y_zp=np.array([0,1]), y_0_inf=1e-4):
    t1=timer()
    x=np.linspace(x_min, x_max, N) #rabi numerov, in za končni izračun
    dx = abs(x[1] - x[0])  # Spacing between x points
    if (n-l)%2==0: y_0_inf=-y_0_inf #sodo št ničel, konvergenca od spodaj
    #iskanje ničle napake za energijo
    try:
        sol = opt.root_scalar(napaka_energija, x0=e0, x1=e0*1.001, method='secant', xtol=10**-10, args=(x_max,l,method,N), maxiter=1000)
        e0 = sol.root #bolj delikatni zač. približki
        print('E_opt:', e0, 'secant', sol.iterations, sol.converged)
        if not sol.converged: raise ValueError('Sekantna metoda ne konvergira.')
    except:
        sol = opt.root(napaka_energija, x0=e0, method='lm', tol=10**-10, args=(x_max,l,method))
        e0 = sol.x[0]
        print('E_opt:', e0, 'nelin', sol.nfev, sol.success)
    t2=timer()
    #final plotting
    if method =='nazaj':
        x=x[::-1]
        y_zp = np.array([y_0_inf, 0]) #zp v neskončnosti
        R, dRdx = solve_ivp(DE, t_span=(x_max,x_min), y0=y_zp, method='RK45', args=(e0,l), t_eval=x).y
        dRdx[-1] = y_zp[1] #lahko težave z odvodi za nazaj v x=0? -> stabilizacija
        dRdx_0 = R[-2]-R[-1]
    elif method == 'naprej':
        R, dRdx = solve_ivp(DE, t_span=(x_min,x_max), y0=y_zp, method='RK45', args=(e0,l), t_eval=x).y
        dRdx_0 = dRdx[0]
    elif method == 'naprej-RK': #konst dx! kot Numerov
        R, dRdx = runge_kutta_4(DE,x,y_zp,e0,l)
        dRdx_0 = dRdx[0]
    elif method == 'numerov':
        R = numerov(x, y_zp, e0, l)
        dRdx=np.ones_like(R)
        dRdx_0 = y_zp[1]
    elif method == 'numerov nazaj':
        x=x[::-1]
        y_zp = np.array([y_0_inf, 0]) #zp v neskončnosti
        R = numerov(x, y_zp, e0, l)
        dRdx=np.ones_like(R)
        dRdx_0 = y_zp[1]
    elif method == 'obojestransko':
        y1, dy1dx = solve_ivp(DE, t_span=(x_min,x[len(x)//2]), y0=y_zp, method='RK45', args=(e0,l), t_eval=x[:len(x)//2+1]).y #naprej
        y_zp = np.array([1e-4, 0]) #zp v neskončnosti
        y2, dy2dx = solve_ivp(DE, t_span=(x_max,x[len(x)//2]), y0=y_zp, method='RK45', args=(e0,l), t_eval=x[:len(x)//2:-1]).y #nazaj
        scaling_factor = y1[-1] / y2[-1]
        y2 *= scaling_factor  # skaliranje, da se točki ujameta
        dy2dx *= scaling_factor 
        R = np.append(y1,y2[::-1]) #združitev
        dRdx = np.append(dy1dx,dy2dx[::-1]) #združitev
        dRdx_0 = dRdx[0]
    if dRdx_0 < 0: # naklon na zač. želimo 1
        R/=-1
        # print('obračam',method)
    norm = np.sqrt(np.sum(np.abs(R)**2 * dx))  # Approximate the integral
    R /= norm  # Normalize the wavefunction
    nzR1 = np.sum((R[:-1]*R[1:]) <0) +1 # ena ničla v izhodišču!
    nzR = count(R)
    if nzR != nzR1: print('ni enako ničel:',nzR,'/=',nzR1)
    if method == 'numerov' or method == 'numerov nazaj':
        return x, R, dRdx, e0, nzR, nzR, t2-t1
    else:
        nzdR1 = np.sum((dRdx[:-1]*dRdx[1:]) <0)
        nzdR = count(dRdx) -1 #ni dodatne ničle
        if nzdR != nzdR1: print('ni enako ničel odvodov:',nzdR,'/=',nzdR1)
        return x, R, dRdx, e0, nzR, nzdR, t2-t1
# n=5 #lupina
# l=1 #podlupina
# N=1000 #št točk x-os
# x_min=1e-9 #singularno v x=0
# y_zp = np.array([0, 1]) #začetna lega, naklon -> naklon NE vpliva na obliko rešitve -> normalizacija
# e0 = -1/n**2 #analitično
# print('E_opt:', e0, 'teoretično')

# method1='naprej-RK' #smer integracije NAPREJ
# x_max = 10*n**1.5 #x_max ocena! naprej: faktor 10, nazaj 20!
# x1,R1,dR1dx,e01,nz1,nzdR1,t = last_f(method1, x_min, x_max, e0, l, N=N)

# method2 = 'nazaj' #smer integracije NAZAJ -> težave za podlupine >0 ?
# x_max = 20*n**1.5 #x_max ocena!
# x2,R2,dR2dx,e02,nz2,nzdR2,t = last_f(method2, x_min, x_max, e0, l, N=N)

# method4 = 'obojestransko' #smer integracije obe smeri do pol
# x_max = 20*n**1.5 #x_max ocena!
# x4,R4,dR4dx,e04,nz4,nzdR4,t = last_f(method4, x_min, x_max, e0, l, N=N)

# method3 = 'numerov'
# x_max = 15*n**1.5 #x_max ocena!
# x3,R3,_,e03,nz3,nzdR3,t = last_f(method3, x_min, x_max, e0, l, N=N)

# method5 = 'numerov nazaj'
# x_max = 20*n**1.5 #x_max ocena!
# x5,R5,_,e05,nz5,nzdR5,t = last_f(method5, x_min, x_max, e0, l, N=N)
# print('pravilno št. prečkanj: ',n-l,', metode:',nz1,nz2,nz4,nz3,nz5) #št. prečkanj x-osi (ničle)
# print('ničle odvodov:',nzdR1,nzdR2,nzdR4,nzdR3,nzdR5)

# x=np.linspace(x_min, max([x1.max(),x2.max(),x3.max()]),N)
# # dx = x[1]-x[0]
# # plt.plot(x, R10(x), label='R10')
# # plt.plot(x, R20(x), label='R20')
# # plt.plot(x, R21(x), label='R21')
# plt.plot(x1, R1, label=f'$e={e01:.1E}$ {method1}')
# plt.plot(x2, R2, '--', label=f'$e={e02:.1E}$ {method2}')
# plt.plot(x4, R4, '--', label=f'$e={e04:.1E}$ {method4}')
# plt.plot(x3, R3, '--', label=f'$e={e03:.1E}$ {method3}')
# plt.plot(x5, R5, '--', label=f'$e={e05:.1E}$ {method5}')
# plt.title(f'$n={n}$, $l={l}$')
# # plt.yscale('log')
# plt.legend()
# plt.show()

# plt.plot(x1, dR1dx, label=f'$e={e01:.1E}$ {method1}')  #odvod = naklon
# plt.plot(x2, dR2dx, '--',label=f'$e={e02:.1E}$ {method2}')
# plt.plot(x4, R4, '--', label=f'$e={e04:.1E}$ {method4}')
# plt.title(f'$n={n}$, $l={l}$')
# plt.legend()
# plt.show()

# prikaz lastnih funkcij
######################################################################
# x_min=1e-9
# dx=0.01

# R=[]
# dRdx=[]
# n_list=[]
# l_list=[]
# e0_list=[]
# x_list=[]
# method = 'obojestransko'
# n_er=0
# # for n in range(1,21): #1. vse
# #     for l in range(0,n):
# for n in [1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50]: #2. le n, l=0
#         l=0
# # for n in [2,3,4,5,6,7]: #3. tirno kvantno
# #     for l in range(1,n):
#         e0=-1/n**2
#         # x_max = 15*n**1.52 #x_max naprej
#         # x_max = 20*n**1.5 #x_max nazaj
#         x_max = 10*n**1.65 #x_max obojestr., numerov
#         N=round((x_max-x_min)/dx)+1
#         xi,Ri,dRidx,e0i,nz,nzdR,t = last_f(method, x_min, x_max, e0, l, N=N)
#         if nz != n-l or nzdR != n-l:
#             n_er+=1
#         if  nzdR != n-l: print('napačne ničle odvoda', nzdR, n-l)
#         n_list.append(n)
#         l_list.append(l)
#         e0_list.append(e0i)
#         R.append(Ri)
#         dRdx.append(dRidx)
#         x_list.append(xi)
# print('izračunov: ',len(e0_list),', napačno ničel: ', n_er)

# fig, ax = plt.subplots(4,4, figsize=(8,6))
# i=j=0
# for k in range(4*4):
#     ax[j,i].plot(x_list[k], R[k], 'k-', label=f'$N={n_list[k]}$, $l={l_list[k]}$', lw=0.75)
#     # ax[j,i].set_xlabel('$x$')
#     ax[j,i].set_ylabel(f'$R_{{{n_list[k]}{l_list[k]}}}(x)$')
#     i+=1
#     if i==4:
#         j+=1
#         i=0
# plt.tight_layout()
# # plt.savefig('./slike/lastne_funk_n.pdf')
# plt.show()

# n_l0=[]
# e0_n=[]
# n1=[]
# l1=[]
# e0_l=[]
# for i in range(len(e0_list)):
#     if l_list[i] == 0:
#         n_l0.append(n_list[i])
#         e0_n.append(e0_list[i])
#     else:
#         n1.append(n_list[i])
#         l1.append(l_list[i])
#         e0_l.append(e0_list[i])
# xx=np.linspace(0.8,22,100)
# plt.plot(xx,1/xx**2,'C0--',label='analitično')
# plt.plot(n_l0, np.abs(e0_n),'ko-',label='$l=0$')
# plt.plot(n1, np.abs(e0_l),'r.',label='$l > 0$')
# plt.yscale('log')
# plt.xlabel('$n$')
# plt.ylabel('$|E|$')
# plt.xticks(np.arange(1,np.max(n_list)+1,2,dtype=int),labels=np.arange(1,np.max(n_list)+1,2,dtype=int))
# plt.legend()
# plt.tight_layout()
# # plt.savefig(f'./slike/energije_n.pdf')
# plt.show()

# odvistnost od dx
######################################################################
# x_min=1e-9
# n=3
# l=2
# R_analit = R32
# e0=-1/n**2

# R_err=[]
# dRdx=[]
# x_list=[]
# dx_list=[]
# energije=[]
# method = 'obojestransko'
# n_er=0
# for dx in [0.5,0.1,0.01,0.001]:
#     # x_max = 15*n**1.52 #x_max naprej
#     # x_max = 20*n**1.5 #x_max nazaj
#     x_max = 10*n**1.65 #x_max obojestr., numerov
#     N=round((x_max-x_min)/dx)+1
#     xi,Ri,dRidx,e0i,nz,nzdR,t = last_f(method, x_min, x_max, e0, l, N=N)
#     if nz != n-l or nzdR != n-l:
#         n_er+=1
#     R_err.append(np.abs(Ri - R_analit(xi)))
#     dRdx.append(dRidx) #TODO error -> odvodi analitične!
#     x_list.append(xi)
#     dx_list.append(dx)
#     energije.append(e0i)
# print('izračunov: ',len(dx_list),', napačno ničel: ', n_er)

# for i in range(len(dx_list)):
#     plt.plot(x_list[i], R_err[i], '-',label=f'{dx_list[i]}')
# plt.yscale('log')
# plt.xlabel('$x$')
# plt.ylabel(f'$|R_{{{n}{l}}} - \mathcal{{R}}_{{{n}{l}}}|$')
# plt.title(f'{method}')
# plt.legend(title='$\\Delta x$')
# # plt.savefig(f'./slike/napaka_dx_n{n}_l{l}_{method}.pdf')
# plt.show()

# plt.plot(dx_list, np.abs(energije),'k.-')
# plt.axhline(1/n**2,c='k',ls='--')
# plt.yscale('log')
# plt.xlabel('$\Delta x$')
# plt.ylabel('$|E|$') #to ne odstopa toliko!!!!
# # plt.savefig(f'./slike/energije_n{n}_l{l}_{method}.pdf')
# plt.show()

# pravilnost rešitev
###################################################################
x_min=1e-9
method = 'obojestransko' #numerov brez @njit -> brez dela linearen ZP
Nn=1+20 #max glavno kvantno št.
dx=0.01

primeri = np.zeros((Nn,Nn))
energije=np.zeros((Nn,Nn))
time=np.zeros((Nn,Nn))
i=n_er=0
for n in range(1,Nn):
    for l in range(0,n):
        try:
            # x_max = 15*n**1.52 #x_max naprej
            # x_max = 20*n**1.5 #x_max nazaj
            x_max = 10*n**1.65 #x_max obojestr., numerov
            N=round((x_max-x_min)/dx)+1
            e0=-1/n**2
            xi,Ri,dRidx,e0i,nz,nzdR,ti = last_f(method, x_min, x_max, e0, l, N=N)
            if nz == n-l and nzdR == n-l:
                primeri[n-1,l] = 1
            else:
                primeri[n-1,l] = -1
                n_er+=1
            i+=1
            energije[n-1,l] = e0i
            time[n-1,l] = ti
        except: pass
print('izračunov: ',i,', napačno ničel: ', n_er)
print(method, 't_avg = ', np.average(time), 't_std = ', np.std(time))

cmap=ListedColormap(['red','white','green'])
bounds = [-1, -0.01, 0.01, 1] 
norm = BoundaryNorm(bounds, cmap.N)
plt.imshow(primeri.T, cmap=cmap, norm=norm)
cbar = plt.colorbar(ticks=[-0.5, 0, 0.5])
cbar.ax.set_yticklabels(['napačno', '', 'pravilno'])
plt.gca().set_yticks(np.arange(0,len(primeri),2,dtype=int))
plt.gca().set_yticklabels(np.arange(0,len(primeri),2,dtype=int))
plt.gca().set_xticks(np.arange(0,len(primeri)-1,2,dtype=int))
plt.gca().set_xticklabels(np.arange(1,len(primeri),2,dtype=int))
plt.grid(False)
plt.xlabel('$n$')
plt.ylabel('$l$')
plt.title(f'$\Delta x ={dx}$, method={method}')
plt.tight_layout()
# plt.savefig(f'./slike/napake_dx{dx}_{method}.pdf')
plt.show()

# povprecen čas simulacij
#dx=0.01
# t_avg1 = [0.022102133333285725, 0.01896191247158933, 0.13799803151921716, 0.14456498095234746, 0.20728056031742584] #old numerov-ZP
# t_std1 = [0.025721475532117653, 0.02417953173690591, 0.17502351237791786, 0.19031919125132438, 0.27049572286197326]
t_avg1 = [0.022102133333285725, 0.01896191247158933, 0.0014167335600286772, 0.0015418532879629394, 0.20728056031742584]
t_std1 = [0.025721475532117653, 0.02417953173690591, 0.0019290054207072038, 0.0020432572494635823, 0.27049572286197326]
metode = ['nazaj', 'naprej', 'numerov', 'numerov nazaj', 'naprej-RK']
#dx=0.1
t_avg2 = [0.022119500226712212, 0.017633315873042735, 0.00029242925176851275, 0.00029448004531173317, 0.024923135600879495]
t_std2 = [00.02584041852950307, 0.022466460159324443, 00.0003325739939212729, 0.0003381553079921072, 0.04405623911170323]
# t_avg2 = [0.022119500226712212, 0.017633315873042735, 0.018161861904804535, 0.018125386394523973, 0.024923135600879495] #old numerov-ZP
# t_std2 = [00.02584041852950307, 0.022466460159324443, 0.023184569603309394, 0.0233436590690293, 0.04405623911170323]
metode = ['nazaj', 'naprej', 'numerov', 'numerov nazaj', 'naprej-RK']

for i in range(len(metode)):
    plt.errorbar(i, 1e3*t_avg2[i], 1e3*t_std2[i], c='b', fmt='.', capsize=5)
    plt.errorbar(i, 1e3*t_avg1[i], 1e3*t_std1[i], c='r', fmt='.', capsize=5, alpha=0.6)
lines=[Line2D([0],[0],c='b',lw=1), Line2D([0],[0],c='r',lw=1)]
plt.xticks(np.arange(0,len(metode),1), labels=metode)
plt.legend(lines, ['0.1', '0.01'], title='$\Delta x$')
plt.ylabel('$t$ [ms]')
plt.yscale('log')
plt.tight_layout()
# plt.savefig(f'./slike/cas.pdf')
plt.show()