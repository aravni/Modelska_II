# -*- coding: utf-8 -*-
# 3. naloga lastne energije schroedingerjeve en.
# 2. del: monokromatska svetloba v snovi - Helmholtzova en.
# ---- povzeto po sprotnem reševanju iz datoteke 3.nal.ipynb ----

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy import optimize as opt
from scipy.integrate import simps
from timeit import default_timer as timer
from numba import njit, jit

plt.style.use('aljaz')

@njit
def lomni(x,Nj=2):
    if x < 1: n = Nj-0.5*x**2
    # elif x < 2: n = 1.25-0.2*x+0.2
    else: n = 1
    return n
    # return np.where(x < 1, 2-0.5*x**2, 1)
@njit
def k_fun(x,λ,k):
    return 1/(4*x**2) + lomni(x)**2*k**2 - λ**2
@njit
def DE(x,y,λ,k):
    R, z = y.reshape((2,len(y)//2))
    dRdx = z
    dzdx = -k_fun(x,λ,k) * R
    return np.append(dRdx,dzdx)
# @njit
def numerov(x, y0, λ, k):
    '''reševanje DE oblike [d²/dx² + k²(x)] y(x) = 0,
    funkcija k(x) naj bo že kvadrirana. Pazi na časovni korak!'''
    h = x[1]-x[0]
    y = np.zeros_like(x)

    y[0] = y0[0] #ZP
    y[1] = y0[0] + h*y0[1] #lin. razvoj

    for i in range(1,len(x)-1):
        a_2 = 1 - h**2*5/12*k_fun(x[i],λ,k)
        a_1 = 1 + h**2*1/12*k_fun(x[i-1],λ,k)
        a_3 = 1 + h**2*1/12*k_fun(x[i+1],λ,k)
        y[i+1] = (2*a_2*y[i] - a_1*y[i-1]) / a_3
    return y
@njit
def runge_kutta_4(f, t, y0, λ, k):
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
        k0 = h*f(ti, y[i], λ, k)
        k1 = h*f(ti + h/2.0, y[i] + k0/2.0, λ, k)
        k2 = h*f(ti + h/2.0, y[i] + k1/2.0, λ, k)
        k3 = h*f(ti + h, y[i] + k2, λ, k)
        y[i+1] = y[i] + (k0 + 2.0*k1 + 2.0*k2 + k3)/6.0
        i+=1
    return np.transpose(y)
@njit
def count(y):
    '''štetje ničel'''
    n=1 #ena ničla na začetku
    for i in range(2,len(y)-1-1):
        if ((y[i]>=0 and y[i+1]<=0 and y[i]!=y[i+1]) or
            (y[i]<=0 and y[i+1]>=0 and y[i]!=y[i+1])):
            n+=1
    return n

# lomni = np.vectorize(lomni)
# x=np.linspace(0,3,1000)
# plt.plot(x,lomni(x))
# plt.xlabel('$x$')
# plt.ylabel('$n(x)$')
# plt.show()

# ročno
#########################################################
k = 1 #prosti parameter (valovno število)
λ = 1.121301781681113 #1.121 #lastna frekvenca
dx=0.1
x_min, x_max = 1e-9, 100
N=round((x_max-x_min)/dx)+1
x = np.linspace(x_min,x_max,N)
y_zp=np.array([0, 0.1]) #ZP (lega, naklon) pri x=0
y_zp_inf=np.array([1e-4, 0]) #ZP (lega, naklon) v neskončnosti

# x=x[::-1] #nazaj bolj stabilno!
# R_rk,dR_rkdx = runge_kutta_4(DE, x, y_zp_inf, λ, k)
# R_n = numerov(x, y_zp_inf, λ, k)

# res = solve_ivp(DE, (x_max,x_min), y0=y_zp_inf, method='RK45', t_eval=x, args=(λ,k)) #nazaj
# R1,dR1dx =res.y #R = lastna oblika
# x1 =res.t
# # res = solve_ivp(DE, (x_min,x_max), y0=y_zp, method='RK45', t_eval=x[::-1], args=(λ,k)) #naprej
# # R2,dR2dx =res.y #R = lastna oblika
# # x2 =res.t

# plt.plot(x1, R1, label='RK45')
# # plt.plot(x2, R2, label='RK45 naprej')
# plt.plot(x[:-1], R_rk[:-1], label='RK4')
# plt.plot(x, R_n, label='Numerov')
# plt.legend()
# plt.show()

# grid search
################################################################################
def napaka_RP(λ, method, k, N=1000, x_max=100, x_min=1e-9):
    y_zp_inf=np.array([1e-4, 0])
    # if k%2 ==0: y_zp_inf[0] *=-1 #sode št ničel!?
    y_zp_0=np.array([0, 0.1])
    if method=='RK45 nazaj':
        R,_ = solve_ivp(DE, (x_max,x_min), y0=y_zp_inf, method='RK45', args=(λ,k)).y
    if method=='RK45 naprej':
        R,_ = solve_ivp(DE, (x_min,x_max), y0=y_zp_0, method='RK45', args=(λ,k)).y
    elif method=='RK4 naprej':
        x=np.linspace(x_min,x_max, N)
        R,_ = runge_kutta_4(DE, x, y_zp_0, λ, k)
    elif method=='RK4 nazaj':
        x=np.linspace(x_max,x_min, N)
        R,_ = runge_kutta_4(DE, x, y_zp_inf, λ, k)
    elif method=='Numerov nazaj':
        x=np.linspace(x_max,x_min, N)
        R = numerov(x, y_zp_inf, λ, k)
    elif method=='Numerov naprej':
        x=np.linspace(x_min,x_max, N)
        R = numerov(x, y_zp_0, λ, k)
    return R[-1] #napaka na RP -> abs
# # lambda_list = np.linspace(0.9, 1.5, 500)
# # lambda_list = np.linspace(3, 6, 1000)
# k=1
# x_max=20
# dx=0.1
# N=round((x_max-x_min)/dx)+1
# lambda_list = np.linspace(k, 2*k, 1000)
# # y_zp_inf=np.array([+1e-4, 0]) #liho
# # y_zp_inf=np.array([-1e-4, 0]) #sodo št. ekstremov -> ni razlike!

# method1='RK45 nazaj'
# t1=timer()
# napaka1=np.zeros_like(lambda_list)
# for i,l in enumerate(lambda_list):
#     napaka1[i] = napaka_RP(l, method1, k, N,x_max)
# t2=timer()
# print(f't_sim = {(t2-t1):.3f} s')
# print(count(napaka1)-1)

# method2='RK4 nazaj'
# t1=timer()
# napaka2=np.zeros_like(lambda_list)
# for i,l in enumerate(lambda_list):
#     napaka2[i] = napaka_RP(l, method2, k, N,x_max)
# t2=timer()
# print(f't_sim = {(t2-t1):.3f} s')
# print(count(napaka2)-1)

# method3='Numerov nazaj'
# t1=timer()
# napaka3=np.zeros_like(lambda_list)
# for i,l in enumerate(lambda_list):
#     napaka3[i] = napaka_RP(l, method3, k, N,x_max)
# t2=timer()
# print(f't_sim = {(t2-t1):.3f} s')
# print(count(napaka3)-1)

# method4='RK45 naprej'
# t1=timer()
# napaka4=np.zeros_like(lambda_list)
# for i,l in enumerate(lambda_list):
#     napaka4[i] = napaka_RP(l, method4, k, N,x_max)
# t2=timer()
# print(f't_sim = {(t2-t1):.3f} s')
# print(count(napaka4)-1)

# method5='RK4 naprej'
# t1=timer()
# napaka5=np.zeros_like(lambda_list)
# for i,l in enumerate(lambda_list):
#     napaka5[i] = napaka_RP(l, method5, k, N,x_max)
# t2=timer()
# print(f't_sim = {(t2-t1):.3f} s')
# print(count(napaka5)-1)

# method6='Numerov naprej'
# t1=timer()
# napaka6=np.zeros_like(lambda_list)
# for i,l in enumerate(lambda_list):
#     napaka6[i] = napaka_RP(l, method6, k, N,x_max)
# t2=timer()
# print(f't_sim = {(t2-t1):.3f} s')
# print(count(napaka6)-1)

# plt.plot(lambda_list, np.abs(napaka1), 'C0-',label=method1)
# plt.plot(lambda_list, np.abs(napaka2), 'C1-',label=method2)
# plt.plot(lambda_list, np.abs(napaka3), 'C2-',label=method3)
# plt.plot(lambda_list, np.abs(napaka4), 'C0--',label=method4)
# plt.plot(lambda_list, np.abs(napaka5), 'C1--',label=method5)
# plt.plot(lambda_list, np.abs(napaka6), 'C2--',label=method6)
# plt.yscale('log')
# plt.legend()
# plt.xlabel('$\lambda$')
# plt.ylabel('abs. napaka RP')
# plt.title(f'$k={k}$')
# # plt.savefig(f'./slike/grid_search_k{k}.pdf')
# plt.show()

# optimalen lambda
#####################################################################################
def bisect(fun, lambda1,lambda2, i,Nmax, roots):
    #poišči interval prehoda ničle -> zgreši če sta na intervalu 2 ničli
    if i<Nmax:
        i+=1
        lambda3=(lambda1+lambda2)/2 #razpolovimo
        f1=fun(lambda1)
        f2=fun(lambda2)
        f3=fun(lambda3)
        if ((f3<0 and f1>0) or (f3>0 and f1<0)):
            bisect(fun,lambda1,lambda3, i,Nmax, roots)
        if ((f3<0 and f2>0) or (f3>0 and f2<0)):
            bisect(fun,lambda3,lambda2, i,Nmax, roots)
    else:
        roots.append([lambda1,lambda2])
def incr(fun, xmin,xmax, N, args):
    '''poišče manjše intervale ničel'''
    roots=[]
    x = np.linspace(xmin,xmax,N)
    for i in range(N-1):
        x1 = x[i]
        x2 = x[i+1]
        f1 = fun(x1, *args)
        f2 = fun(x2, *args)
        if (f1<0 and f2>0) or (f1>0 and f2<0): roots.append([x1,x2])
    return roots
def nicle(method, k, N=100, Ndiv=100):
    #iskanje manjših intervalov ničel, vedno med k in 2*k
    roots=incr(napaka_RP, k, 2*k, N=Ndiv, args=(method, k, N, x_max))

    roots_1=[] #točnejša ničla na vsakem intervalčku
    for root in roots:
        try:
            # sol = opt.root_scalar(napaka_RP, x0=root[0], x1=root[1], method='secant', xtol=10**-15, args=(method,k,N,x_max), maxiter=1000)
            sol = opt.root_scalar(napaka_RP, bracket=(root[0],root[1]), method='brentq', xtol=1e-15, args=(method, k, N, x_max), maxiter=1000)
            lambda0 = sol.root #bolj delikatni zač. približki
            print('lambda_opt:', lambda0, 'brentq', sol.iterations, sol.converged)
            if not sol.converged: raise ValueError('Brent metoda ne konvergira.')
        except ValueError:
            # sol = opt.root(napaka_RP, x0=root[0], method='broyden1', tol=10**-5, args=(method,k,N,x_max), options={'maxiter':1000})
            sol = opt.root(napaka_RP, x0=root[0], method='hybr', tol=10**-5, args=(method,k,N,x_max), options={'maxiter':1000})
            # sol = opt.root(napaka_RP, x0=root[0], method='lm', tol=10**-10, args=(method,k,N,x_max))
            lambda0 = sol.x[0]
            print('lambda_opt:', lambda0, 'nelin', sol.nfev, sol.success)
        roots_1.append(lambda0)
    return roots_1
def last_f(method,lambda0,k,x_min,x_max,N):
    y_zp_inf=np.array([1e-4,0])
    y_zp_0=np.array([0, 0.1])
    x,dx=np.linspace(x_max,x_min, N, retstep=True)
    if method=='RK45 nazaj':
        R,_ = solve_ivp(DE, (x_max,x_min), y0=y_zp_inf, method='RK45', t_eval=x, args=(lambda0,k)).y
        R=R[::-1]
    elif method=='RK45 naprej':
        R,_ = solve_ivp(DE, (x_min,x_max), y0=y_zp_0, method='RK45', t_eval=x[::-1], args=(lambda0,k)).y
    elif method=='RK4 nazaj':
        R,_ = runge_kutta_4(DE, x, y_zp_inf, lambda0, k)
        R=R[::-1] #nestabilno blizu x=0
    elif method=='RK4 naprej':
        R,_ = runge_kutta_4(DE, x[::-1], y_zp_0, lambda0, k)
    elif method=='Numerov nazaj':
        R = numerov(x, y_zp_inf, lambda0, k)[::-1]
    elif method=='Numerov naprej':
        R = numerov(x[::-1], y_zp_0, lambda0, k)
    #normiranje
    norm = np.sqrt(np.sum(np.abs(R)**2 * abs(dx)))
    R /= norm
    if R[1] < R[0]: R *= -1 #obrnemo
    nz=count(R) #št. ničel  = št.ekstremov
    return R, nz

k=1
dx=0.01
method='Numerov nazaj'
x_max=10/(0.5*k)
Ndiv=10 #št delitev za ničle
N=round((x_max-x_min)/dx)+1
x=np.linspace(x_min,x_max,N)

# roots = nicle(method, k, N, Ndiv)
# for i in range(len(roots)):
#     R, nz = last_f(method, roots[i],k, x_min,x_max,N)
#     plt.plot(x,R,label=f'$\lambda={roots[i]:2.3f}$; $N={nz:n}$')
# plt.legend(title=f'$\lambda$; št.vezanih stanj: {len(roots)}')
# plt.xlabel('$x$')
# plt.ylabel('$R(x)$')
# plt.title(f'{method}')
# # plt.savefig(f'./slike/modes_k{k}_{method}.pdf')
# plt.show()

# disperzijska relacija
#################################################################################
# dx=0.01
# Ndiv=100
# method='Numerov nazaj'

# k_values = np.linspace(0.8, 10, 50)
# lambda_vals = {}
# time=[]
# nz_list=[]
# k_dict={}
# for k in k_values:
#     print('k = ',k)
#     x_max = 20
#     N=round((x_max-x_min)/dx)+1
#     t1=timer()
#     roots = nicle(method, k, N, Ndiv)
#     t1=timer()-t1
#     time.append(t1)
#     for i in range(len(roots)):
#         Ri,nzi = last_f(method,roots[i],k,x_min,x_max,N)
#         if nzi not in lambda_vals: #novo mode valovanja (št. ekstremov)
#             lambda_vals[nzi] = []
#             k_dict[nzi] = []
#             nz_list.append(nzi)
#         lambda_vals[nzi].append(roots[i])
#         k_dict[nzi].append(k)
# nz_list=np.sort(nz_list)

# for i in range(len(nz_list)):
#     plt.plot(k_dict[nz_list[i]], lambda_vals[nz_list[i]], '.-', label=f'{nz_list[i]:n}')
# xx=np.linspace(0,np.max(k_values),100)
# plt.plot(xx,xx,'k--',label='$\lambda=k$')
# plt.xlabel('$k$')
# plt.ylabel('$\lambda$')
# plt.legend(title='vezana stanja')
# # plt.savefig(f'./slike/disperzijska_2.pdf')
# plt.show()

# plt.plot(k_values, time, 'k.-')
# plt.xlabel('$k$')
# plt.ylabel('$t$ [s]')
# # plt.savefig('./slike/cas_helmholtz.pdf')
# plt.show()

############### single mode
# s_mode=min(k_dict[2])
# print(f'single mode fiber do k = {s_mode}')
Nj = [1.5, 1.75, 2, 2.25, 2.5, 2.75, 3]
single_mode = [5.118367346938775, 3.6163265306122447, 2.8653061224489793, 2.3020408163265307, 2.114285714285714, 1.9265306122448977, 1.7387755102040816]

plt.plot(Nj, single_mode, 'k.-')
plt.xlabel('$N_j$')
plt.ylabel('$k_{single-mode}$')
plt.show()

############### presečišče
# coef=[]
# for i in range(len(nz_list)):
#     x = k_dict[nz_list[i]]
#     y = lambda_vals[nz_list[i]]
#     k,n = np.polyfit(x, y, 1)
#     coef.append((k,n))

# A = []
# B = []
# for (k,n) in coef: # Rearranging the line equation to m*x - y = -b
#     A.append([k,-1])
#     B.append(-n)
# A = np.array(A)
# B = np.array(B)
# x_i, y_i = np.linalg.lstsq(A, B, rcond=None)[0]
# print(f'Presečišče: ({x_i}, {y_i})')

presecisca = [[1.5,-8.27,-14.023], [1.75,48.305,83.217], [2,24.775,48.671], [2.25,20.245,45.052], [2.5,19.93,51.55], [2.75,20.288,61.488],[3,6.84,10.91]] #odvisno od N v n(x,N)!! navodila: 2; [N,x,y]
presecisca=np.array(presecisca)
for i in range(presecisca.shape[0]):
    plt.plot(presecisca[i,1], presecisca[i,2], 'x', label=f'{presecisca[i,0]}')
plt.legend(title='$N_j$')
plt.xlabel('$k$')
plt.ylabel('$\lambda$')
# plt.savefig(f'./slike/presecisce.pdf')
plt.show()