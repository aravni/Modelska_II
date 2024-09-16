# -*- coding: utf-8 -*-
# nihalo z N diskretnimi masami
# algoritem aproksimira neraztegljivo vrvico z Euler-Lagrangeevimi gibalnimi en., predavanja Modelska II pa 2.Newtonov zakon
# analiza periodike nihanja

import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer
from scipy.integrate import solve_ivp
import matplotlib.animation as animation
from scipy.signal import windows, find_peaks
from alignYaxes import *
from ProcessPriority import *
from numba import jit
from matplotlib.ticker import MaxNLocator
from concurrent.futures  import ProcessPoolExecutor
from scipy.optimize import curve_fit

HighPriority()
plt.style.use('aljaz')
np.printoptions(linewidth=np.nan)

@jit(nopython=True)
def post_process(phi, phi_dot, m, L, g, Nt):
    N=len(m)
    l=L/N
    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)
    x = np.zeros((Nt, N))
    y = np.zeros((Nt, N))
    Ek = np.zeros(Nt)
    Ep = np.zeros(Nt)
    for t_i in range(Nt):
        x[t_i, :] = np.cumsum(sin_phi[t_i, :] * l)
        y[t_i, :] = -np.cumsum(cos_phi[t_i, :] * l)
        Ej = np.zeros((2, N))
        for i in range(N):
            Ej[0, i] = np.sum(phi_dot[t_i, :i+1] * cos_phi[t_i, :i+1])
            Ej[1, i] = np.sum(phi_dot[t_i, :i+1] * sin_phi[t_i, :i+1])
            Ep[t_i] += m[i] * g * (L - l * cos_phi[t_i, i])
            Ek[t_i] += 0.5 * m[i] * l**2 * (Ej[0, i]**2 + Ej[1, i]**2)
    E = Ek + Ep
    Ek /= E[0]
    Ep /= E[0]
    E /= E[0]
    return x,y, Ek,Ep,E

N = 10 #št. diskretizacije
L = 1. #m, dolžina vrvice
M = 1. #kg, masa vrvice
d = 0.007 #dušenje razmernik

g = 1. #m/s²
l = L/N #konstantna diskretizacija!

m = np.zeros(N)
m[:] = M/N # 1. konst. masa po vrvi
# m[:] = 1e-3 #2. aprox. mat. nihalo
# m[-1] = M*100000

ω_0 = np.sqrt(g/L) #nedušeno mat. nihalo
ω_0d = ω_0 * np.sqrt(1-d**2) #dušeno mat. nihalo
T = 50 * 2*np.pi/ω_0d #fiz čas - teoretično št. nihajev, s
dt = 0.01 #s
t = np.arange(0,T,dt)

#gibalne enačbe
def GE(t, y):
    phi, phi_dot = y.reshape((2,N))
    A = np.zeros((N,N))
    b = np.zeros((N))
    for i in range(N):
        bi = 0
        for j in range(N):
            c = N - max(i,j) + 1
            A[i,j] = c*np.cos(phi[i]-phi[j]) #bolj hitro reševanje, če sistem delimo z 'l^2'
            bi += c*np.sin(phi[i]-phi[j])*phi_dot[j]**2
        d_c = 2*A[i,i]*ω_0/l #kritično dušenje
        b[i] = -bi - g/l*(N-i+1)*np.sin(phi[i]) - d*d_c*phi_dot[i] #viskozno dušenje
    return np.concatenate((phi_dot, np.linalg.solve(A,b)))

if __name__ == '__main__':
    #################################
    #ZP
    phi_0 = np.ones((N)) * np.deg2rad(45)
    phi_dot_0 = np.zeros((N)) #rad/s
    y0 = np.concatenate((phi_0, phi_dot_0))

    t1 = timer()
    res = solve_ivp(GE,(0,T),y0,method='RK45',t_eval=t)
    t2 = timer()
    print(f't_sim = {(t2-t1):.1f}')
    y_res = res.y.T
    phi = y_res[:,:N]
    phi_dot = y_res[:,N:]

    # POST koordinate, energije
    # x=[]
    # y=[]
    # Ek= []
    # Ep=[]
    # for t_i in range(len(t)):
    #     xi = np.zeros((N))
    #     yi = np.zeros((N))
    #     Ei=Epi=0
    #     for i in range(N):
    #         xi[i] = xi[i-1] + np.sin(phi[t_i,i])*l
    #         yi[i] = yi[i-1] + np.cos(phi[t_i,i])*l
    #         Ej=[0,0]
    #         for j in range(i+1): #2. osnovna def.
    #             Ej[0] += phi_dot[t_i,j]*np.cos(phi[t_i,j])
    #             Ej[1] += phi_dot[t_i,j]*np.sin(phi[t_i,j])
    #             Epi += m[i]*g*(L-l*np.cos(phi[t_i,j]))
    #         Ei += 0.5*m[i]*l**2*(Ej[0]**2 + Ej[1]**2)            
    #     x.append(xi)
    #     y.append(-yi)
    #     Ek.append(Ei) #translacija
    #     Ep.append(Epi)
    # x = np.array(x)
    # y = np.array(y)
    # E = np.array(Ek)+np.array(Ep)
    # Ek=np.array(Ek)/E[0]
    # Ep=np.array(Ep)/E[0]
    # E=E/E[0] #normiranje

    # FASTER POST LOOP !!!
    x,y, Ek,Ep,E = post_process(phi, phi_dot, m,L,g,len(t))
    print('t_post = ', timer()-t2)

    # FFT
    ##################################
    Nt = len(t)
    dist = int(ω_0d/2/np.pi * (1/dt))
    okno = windows.hann(Nt)
    f = np.fft.rfftfreq(Nt, d=dt)

    # ############### 1. prva masa
    # index=0
    # X = np.fft.rfft((x[:,index]-np.average(x[:,index]))*okno)*2/Nt
    # Y = np.fft.rfft((y[:,index]-np.average(y[:,index]))*okno)*2/Nt
    # Ep_f = np.fft.rfft(Ep*okno)*2/Nt
    # Ek_f = np.fft.rfft(Ek*okno)*2/Nt
    # E_f = np.fft.rfft(E*okno)*2/Nt
    # peaks_x, dict = find_peaks(np.abs(X), height=1e-6, distance=dist, prominence=1e-4)
    # heights_x = dict['peak_heights']
    # peaks_y, dict = find_peaks(np.abs(Y), height=1e-6, distance=dist, prominence=1e-4)
    # heights_y = dict['peak_heights']
    # peaks_Ep, dict = find_peaks(np.abs(Ep_f), height=1e-6, distance=dist, prominence=1e-4)
    # heights_Ep = dict['peak_heights']
    # peaks_Ek, dict = find_peaks(np.abs(Ek_f), height=1e-6, distance=dist, prominence=1e-4)
    # heights_Ek = dict['peak_heights']
    # peaks_E, dict = find_peaks(np.abs(E_f), height=1e-6, distance=dist, prominence=1e-4)
    # heights_E = dict['peak_heights']

    # try: N_p=int(max(peaks_x.max(), peaks_y.max())*2)
    # except: N_p=Nt//1000
    # plt.axvline(ω_0d/2/np.pi, c='k', ls='--', label='$\\omega_{0d}$')
    # plt.axvline(ω_0d/np.pi, c='k', ls=':', label='$2\cdot\\omega_{0d}$')
    # plt.plot(f[:N_p],np.abs(X)[:N_p], 'C0',label=f'$x$')
    # plt.plot(f[:N_p],np.abs(Y)[:N_p], 'C2', label=f'$y$')
    # plt.plot(f[:N_p],np.abs(Ep_f)[:N_p], 'C1', label=f'$E_p$')
    # plt.plot(f[:N_p],np.abs(Ek_f)[:N_p], 'C3', label=f'$E_k$')
    # plt.plot(f[:N_p],np.abs(E_f)[:N_p], 'C4--', label=f'$E$')
    # plt.plot(f[peaks_x], heights_x, 'C0o')
    # plt.plot(f[peaks_y], heights_y, 'C2o')
    # plt.plot(f[peaks_Ep], heights_Ep, 'C1o')
    # plt.plot(f[peaks_Ek], heights_Ek, 'C3o')
    # plt.plot(f[peaks_E], heights_E, 'C4o')

    # plt.legend()
    # plt.yscale('log')
    # plt.xlabel('$f$ [Hz]')
    # plt.ylabel('$|X|$')
    # plt.title(f'$\\varphi_0 = {np.rad2deg(phi_0[0]):n}, \\zeta={d:.0E}, N_{{all}}={N}, N_i={index+1}$')
    # plt.tight_layout()
    # plt.show()

    # ############### 2. zadnja masa
    # index=N-1
    # X = np.fft.rfft((x[:,index]-np.average(x[:,index]))*okno)*2/Nt
    # Y = np.fft.rfft((y[:,index]-np.average(y[:,index]))*okno)*2/Nt
    # Ep_f = np.fft.rfft(Ep*okno)*2/Nt
    # Ek_f = np.fft.rfft(Ek*okno)*2/Nt
    # E_f = np.fft.rfft(E*okno)*2/Nt
    # peaks_x, dict = find_peaks(np.abs(X), height=1e-6, distance=dist, prominence=1e-6)
    # heights_x = dict['peak_heights']
    # peaks_y, dict = find_peaks(np.abs(Y), height=1e-6, distance=dist, prominence=1e-6)
    # heights_y = dict['peak_heights']
    # peaks_Ep, dict = find_peaks(np.abs(Ep_f), height=1e-6, distance=dist, prominence=1e-6)
    # heights_Ep = dict['peak_heights']
    # peaks_Ek, dict = find_peaks(np.abs(Ek_f), height=1e-6, distance=dist, prominence=1e-6)
    # heights_Ek = dict['peak_heights']
    # peaks_E, dict = find_peaks(np.abs(E_f), height=1e-6, distance=dist, prominence=1e-6)
    # heights_E = dict['peak_heights']

    # try: N_p=int(max(peaks_x.max(), peaks_y.max())*2)
    # except: N_p=Nt//1000
    # plt.axvline(ω_0d/2/np.pi, c='k', ls='--', label='$\\omega_{0d}$')
    # plt.axvline(ω_0d/np.pi, c='k', ls=':', label='$2\cdot\\omega_{0d}$')
    # plt.plot(f[:N_p],np.abs(X)[:N_p], 'C0',label=f'$x$')
    # plt.plot(f[:N_p],np.abs(Y)[:N_p], 'C2', label=f'$y$')
    # plt.plot(f[:N_p],np.abs(Ep_f)[:N_p], 'C1', label=f'$E_p$')
    # plt.plot(f[:N_p],np.abs(Ek_f)[:N_p], 'C3', label=f'$E_k$')
    # plt.plot(f[:N_p],np.abs(E_f)[:N_p], 'C4--', label=f'$E$')
    # plt.plot(f[peaks_x], heights_x, 'C0o')
    # plt.plot(f[peaks_y], heights_y, 'C2o')
    # plt.plot(f[peaks_Ep], heights_Ep, 'C1o')
    # plt.plot(f[peaks_Ek], heights_Ek, 'C3o')
    # plt.plot(f[peaks_E], heights_E, 'C4o')

    # plt.legend()
    # plt.yscale('log')
    # plt.xlabel('$f$ [Hz]')
    # plt.ylabel('$|X|$')
    # plt.title(f'$\\varphi_0 = {np.rad2deg(phi_0[0]):n}, \\zeta={d:.0E}, N_{{all}}={N}, N_i={index+1}$')
    # plt.tight_layout()
    # plt.show()

    # ################## 3. lastne freq po dolžini
    # ang_x=[]
    # ang_y=[]
    # px=[]
    # py=[]
    # peaks_x=np.zeros((N,5)) #max 5 lastnih frek
    # peaks_y=np.zeros((N,5)) #max 5 lastnih frek
    # for i in range(N):
    #     X=np.fft.rfft((x[:,i]-np.average(x[:,i]))*okno)*2/Nt
    #     Y=np.fft.rfft((y[:,i]-np.average(y[:,i]))*okno)*2/Nt
    #     ang_x.append(np.angle(X))
    #     ang_y.append(np.angle(Y))
    #     X = np.abs(X)
    #     Y = np.abs(Y)
    #     max_px = np.argmax(X[5:])
    #     max_py = np.argmax(Y[5:])
    #     px.append(f[5+max_px]) #le prva freq
    #     py.append(f[5+max_py])
    #     peaks, _ = find_peaks(X, height=1e-6, distance=dist, prominence=1e-6)
    #     p_xi = peaks[peaks>5]
    #     peaks, _ = find_peaks(Y, height=1e-6, distance=dist, prominence=1e-6)
    #     p_yi = peaks[peaks>5]
    #     Nf = min([len(p_xi),5])
    #     peaks_x[i,:Nf] = f[p_xi[:Nf]] #prvih Nf lastnih f
    #     Nf = min([len(p_yi),5])
    #     peaks_y[i,:Nf] = f[p_yi[:Nf]] #prvih Nf lastnih f

    # # fig,ax=plt.subplots()
    # # ax.plot(px, 'C0.-', label='$x$')
    # # # ax.axhline(ω_0d/2/np.pi, c='C0', ls='--', label='$f_{0d}$')
    # # ax.set_ylabel('x-kordinata $f$', c='C0')
    # # ax.set_xlabel('$N_i$')
    # # ax.tick_params(axis='y', labelcolor='C0')
    # # ax1=plt.twinx(ax)
    # # ax1.plot(py, 'C1.-', label='$y$')
    # # # ax1.axhline(ω_0d/np.pi, c='C1', ls='--', label='$2\cdot f_{0d}$')
    # # ax1.set_ylabel('y-kordinata $f$', c='C1')
    # # ax1.tick_params(axis='y', labelcolor='C1')
    # # ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    # # alignYaxes([ax,ax1])
    # # plt.title(f'$\\varphi_0 = {np.rad2deg(phi_0[0]):n}, \\zeta={d:.0E}, N_{{all}}={N}$')
    # # plt.tight_layout()
    # # plt.show()

    # fig,ax=plt.subplots()
    # c=['C0', 'C1', 'C2', 'C3', 'C4']
    # for i in range(peaks_y.shape[1]):
    #     if np.any(peaks_x[:, i] != 0) and np.any(peaks_y[:, i] != 0):  # Only plot if not all values are zero
    #         plt.plot(peaks_x[:,i],'.-', c=c[i], label=f'{i+1}.')
    #         plt.plot(peaks_y[:,i],'.--', c=c[i])
    # ax.axhline(ω_0d/2/np.pi, c='k', ls=':', label='$f_{0d}$')

    # plt.legend(title='lastne frekvence')
    # plt.title('polna črta: x-smer, prekinjena: y-smer')
    # plt.xlabel('$N_i$')
    # plt.ylabel('$f$ [Hz]')
    # ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    # plt.show()

    # # fazni zamik v točki
    # i=N-1
    # plt.plot(f[:1000], ang_x[i][:1000], '-', label='$x$', alpha=0.6)
    # plt.plot(f[:1000], ang_y[i][:1000], '-', label='$y$', alpha=0.6)
    # plt.axvline(ω_0d/2/np.pi, c='k', ls=':', label='$f_{0d}$')
    # plt.axvline(ω_0d/np.pi, c='k', ls=':', label='$2\cdot f_{0d}$')
    # plt.yticks([-np.pi,0,np.pi], labels=['$-\pi$', '$0$', '$\pi$'])
    # plt.xlabel('$f$ Hz')
    # plt.ylabel('fazni zamik $\\varphi$')
    # plt.legend(loc='upper right')
    # plt.tight_layout()
    # plt.show()

    # prva lastna freq odvisno od diskretizacije on ZP
    ###########################################################

    N_list = np.linspace(1, 20, 10, dtype=int)
    x_list=np.zeros((len(N_list)))
    print('simulacija za diskretizacije N:')
    input(N_list)
    for i,N in enumerate(N_list):
        t1 = timer()
        l = L/N #konstantna diskretizacija!
        m = np.zeros(N)
        m[:] = M/N # 1. konst. masa po vrvi

        phi_0 = np.ones((N)) * np.deg2rad(45)
        phi_dot_0 = np.zeros((N)) #rad/s
        y0 = np.concatenate((phi_0, phi_dot_0))

        res = solve_ivp(GE,(0,T),y0,method='RK45',t_eval=t)
        y_res = res.y.T
        phi = y_res[:,:N]
        phi_dot = y_res[:,N:]

        x,y, Ek,Ep,E = post_process(phi, phi_dot, m,L,g,len(t))
        X = np.abs(np.fft.rfft((x[:,-1]-np.average(x[:,-1]))*okno)*2/Nt) #zadnja masa!
        x_list[i] = f[5+np.argmax(X[5:])] #1. lastna
        t2 = timer()
        print(f'N = {N}, t_sim = {(t2-t1):.1f}')

    def fun(x, A, B):
        return A + B* np.log10(x) #np.emath.logn(base, x)
    def fun(x, A, B, C): #sigmoid
        return A / (1 + B * np.exp(-C * x))
    (A, B, C),_ = curve_fit(fun, N_list, x_list, p0=(1,1,1), method='lm', maxfev=10000)
    xx=np.inf
    x_ex = fun(xx, A, B,C)
    x=np.linspace(N_list[0], N_list[-1],100)
    x_f=fun(x,A,B,C)
    plt.plot(N_list, x_list, '.-')
    plt.axhline(ω_0d/2/np.pi, c='k', ls=':', label='$f_{0d}$')
    plt.plot(x,x_f, 'k-', label=f'$f(x) = \\frac{{{A}}}{{1+{B} \exp(-{C} x)}}$')
    plt.axhline(x_ex, c='k', ls='--', label=f'$f(\infty)$')
    plt.xlabel('$N_{all}$')
    plt.ylabel('1. lastna $f$ [Hz]')
    plt.legend(loc='upper left', bbox_to_anchor=(0,0.9))
    plt.tight_layout()
    plt.xscale('log')
    plt.show()

    # zač kot v deg
    N=10
    l = L/N #konstantna diskretizacija!
    m = np.zeros(N)
    m[:] = M/N # 1. konst. masa po vrvi
    phi_list = np.linspace(0,180, 10, dtype=int)
    x_list=np.zeros((len(phi_list),N))
    print('simulacija za začetne kote v deg:')
    input(phi_list)
    for i,phi_i in enumerate(phi_list):
        t1 = timer()
        phi_0 = np.ones((N)) * np.deg2rad(phi_i)
        phi_dot_0 = np.zeros((N)) #rad/s
        y0 = np.concatenate((phi_0, phi_dot_0))

        res = solve_ivp(GE,(0,T),y0,method='RK45',t_eval=t)
        y_res = res.y.T
        phi = y_res[:,:N]
        phi_dot = y_res[:,N:]

        x,y, Ek,Ep,E = post_process(phi, phi_dot, m,L,g,len(t))
        for j in range(N):
            X = np.abs(np.fft.rfft((x[:,j]-np.average(x[:,j]))*okno)*2/Nt)
            x_list[i,j] = f[5+np.argmax(X[5:])] #1. lastna
        t2 = timer()
        print(f'phi_0 = {phi_i}, t_sim = {(t2-t1):.1f}')

    for i in [0,N//2,N]:
        plt.plot(phi_list, x_list[:,i], '.-', label=f'$N_i={i+1}$')
    plt.axhline(ω_0d/2/np.pi, c='k', ls=':', label='$f_{0d}$')
    plt.legend(title='delec na vrvi')
    plt.xlabel('$\\varphi_0$ [deg]')
    plt.ylabel('1. lastna $f$ [Hz]')
    plt.tight_layout()
    plt.show()
