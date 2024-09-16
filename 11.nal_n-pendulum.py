# -*- coding: utf-8 -*-
# nihalo z N diskretnimi masami
# algoritem aproksimira neraztegljivo vrvico z Euler-Lagrangeevimi gibalnimi en., predavanja Modelska II pa 2.Newtonov zakon

import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer
from scipy.integrate import solve_ivp
import matplotlib.animation as animation
from numba import jit, njit
from ProcessPriority import *

HighPriority()
plt.style.use('aljaz')
np.printoptions(linewidth=np.nan)

@njit
def GE_compile(t, y, N, d, g, l, omega_0):
    '''gibalne enačbe'''
    phi, phi_dot = y.reshape((2,N))
    A = np.zeros((N,N))
    b = np.zeros((N))
    for i in range(N):
        bi = 0
        for j in range(N):
            c = N - max(i,j) + 1
            A[i,j] = c*np.cos(phi[i]-phi[j]) #bolj hitro reševanje, če sistem delimo z 'l^2'
            bi += c*np.sin(phi[i]-phi[j])*phi_dot[j]**2
        d_c = 2*A[i,i]*omega_0/l #kritično dušenje
        b[i] = -bi - g/l*(N-i+1)*np.sin(phi[i]) - d*d_c*phi_dot[i] #viskozno dušenje
    return np.concatenate((phi_dot, np.linalg.solve(A,b)))
@njit
def RK4(f, t, y, h, N, d, g, l, omega_0):
    '''Perform a single RK4 step.'''
    k0 = h*f(t, y, N, d, g, l, omega_0)
    k1 = h*f(t + h/2.0, y + k0/2.0, N, d, g, l, omega_0)
    k2 = h*f(t + h/2.0, y + k1/2.0, N, d, g, l, omega_0)
    k3 = h*f(t + h, y + k2, N, d, g, l, omega_0)
    return (k0 + 2.0*k1 + 2.0*k2 + k3)/6.0
@njit
def runge_kutta_4(f, t, y0, d, g, L, omega_0):
    """
    Metoda Runge-Kutta 4. reda za reševanje diferencialne enačbe: y' = f(t, y)
    
    :param f: funkcija, ki jo kličemo s parametroma t in y in vrne 
              vrednost prvega odvoda
    :param t: ekvidistanten časovni vektor oz. neodvisna spremenljivka
    :param y0: začetna vrednost
    """    
    y = np.zeros(shape=(len(t),len(y0)))
    y[0] = y0
    h = t[1]-t[0]
    N=len(y0)//2
    l=L/N

    for i in range(len(t)-1):
        y[i+1] = y[i] + RK4(f, t[i], y[i], h, N, d, g, l, omega_0)

    return y
@njit
def post_process(phi, phi_dot):
    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)
    x = np.zeros((len(t), N))
    y = np.zeros((len(t), N))
    Ek = np.zeros(len(t))
    Ep = np.zeros(len(t))
    for t_i in range(len(t)):
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

hist = 700 #zgodovina za prikaz
g = 1. #m/s²
l = L/N #konstantna diskretizacija!

m = np.zeros(N)
m[:] = M/N # 1. konst. masa po vrvi
# m[:] = 1e-3 #2. aprox. mat. nihalo
# m[-1] = M*100000

ω_0 = np.sqrt(g/L) #nedušeno mat. nihalo
ω_0d = ω_0 * np.sqrt(1-d**2) #dušeno mat. nihalo
T = 5 * 2*np.pi/ω_0d #fiz čas - teoretično št. nihajev, s
dt = 0.01 #s; za animacijo je OK 0.01
t = np.arange(0,T,dt)

#gibalne enačbe
# @njit
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
    # return GE_compile(t,y,N,d,g,l,ω_0)
def GE1(t, y): #NI OK!! po Yesilyurt
    phi, phi_dot = y.reshape((2,N))
    A = np.zeros((N,N))
    b = np.zeros((N))
    for i in range(N):
        for j in range(N):
            mk=np.sum(m[j:])
            if j!=i: A[i,j] += mk*l*np.cos(phi[i]-phi[j]) #ok
            b[i] -= mk*l*np.sin(phi[i]-phi[j])*phi_dot[i]*phi_dot[j] #nelin! ok
            b[i] -= mk*l*np.sin(phi[j]-phi[i])*(phi_dot[i]-phi_dot[j])*phi_dot[j] #nelin! ok
            if j>=i:
                b[i] -= m[j]*g*np.sin(phi[i]) #ok
                A[i,i] += m[j]*l #ok
    try: phi_ddot = np.linalg.solve(A, b)  # Solve for angular accelerations
    except np.linalg.LinAlgError:
        A += np.eye(N) * 1e-5
        phi_ddot = np.linalg.solve(A, b)
    return np.concatenate((phi_dot, phi_ddot))

#################################
#ZP - 1
phi_01=phi_0 = np.ones((N)) *np.deg2rad(-180) #konst. naklon, rad
phi_dot_0 = np.zeros((N)) #rad/s
y0 = np.concatenate((phi_0, phi_dot_0))

t1 = timer()
res = solve_ivp(GE,(0,T),y0,method='RK45',t_eval=t)
t2 = timer()
print(f't_sim = {(t2-t1):.1f}')

y_res = res.y.T
phi = y_res[:,:N]
phi_dot = y_res[:,N:]

# t1 = timer() #najhitreje scipy brez numba, drugače hitrejša lastna koda z numba!!
# r = runge_kutta_4(GE_compile, t, y0, d, g,L,ω_0)
# t2 = timer()
# print(f't_sim = {(t2-t1):.1f}')

# POST koordinate, energije
# x=[]
# y=[]
# Ekt=[]
# Ekr=[]
# Ep=[]
# for t_i in range(len(t)):
#     xi = np.zeros((N))
#     yi = np.zeros((N))
#     Ei=Epi=0
#     for i in range(0,N):
#         xi[i] = xi[i-1] + np.sin(phi[t_i,i])*l
#         yi[i] = yi[i-1] + np.cos(phi[t_i,i])*l
#         # for j in range(N): #1. končna definicija = sist. enačb
#         #     Ei += m[i]*0.5*phi_dot[t_i,j]**2*l**2
#         #     if i!=j: Ei += m[i]*l**2*phi_dot[t_i,i]*phi_dot[t_i,j]*np.cos(phi_dot[t_i,i]-phi_dot[t_i,j])
#         # Epi += g*m[i]*(L-(N-i+1)*np.cos(phi[t_i,i])*l)
#         Ej=[0,0]
#         for j in range(i+1): #2. osnovna def.
#             Ej[0] += phi_dot[t_i,j]*np.cos(phi[t_i,j])
#             Ej[1] += phi_dot[t_i,j]*np.sin(phi[t_i,j])
#             Epi += m[i]*g*(L-l*np.cos(phi[t_i,j]))
#         Ei += 0.5*m[i]*l**2*(Ej[0]**2 + Ej[1]**2)           
#         # for j in range(N): #3. Yesilyurt
#         #     for k in range(N):
#         #         Ei += 0.5*m[i]*l*l*phi_dot[t_i,j]*phi_dot[t_i,k]*np.cos(phi[t_i,j]-phi[t_i,k])
#     x.append(xi)
#     y.append(-yi)
#     Ekt.append(Ei) #translacija
#     # Ep.append(g*np.sum(m*(L-yi)))
#     Ep.append(Epi)
# E1 = np.array(Ekt)+np.array(Ep)
# x1 = np.array(x)
# y1 = np.array(y)
# Ekt1= np.array(Ekt)/E1[0]
# Ep1= np.array(Ep)/E1[0]
# E1=E1/E1[0] #normiranje
x1,y1, Ekt1,Ep1,E1 = post_process(phi, phi_dot)
print('t_post = ', timer()-t2)

#################################
# m[:] = M/N # 1. konst. masa po vrvi
# m[:] = 1e-3 #2. aprox. mat. nihalo
# m[-1] = M*100000

#ZP - 2
phi_02=phi_0 = np.ones((N)) *np.deg2rad(10) #konst. naklon, rad
# phi_02=phi_0 = l*np.linspace(0,N-1,N) +np.deg2rad(45) #parabola, rad
phi_dot_0 = np.zeros((N)) #rad/s
y0 = np.concatenate((phi_0, phi_dot_0))

t1 = timer()
res = solve_ivp(GE,(0,T),y0,method='RK45',t_eval=t)
# y_res = runge_kutta_4(GE, t, y0)
t2 = timer()
print(f't_sim = {(t2-t1):.1f}')

y_res = res.y.T
phi = y_res[:,:N]
phi_dot = y_res[:,N:]

# POST koordinate, energije
x=[]
y=[]
Ekr=[]
Ekt= []
Ep=[]
for t_i in range(len(t)):
    xi = np.zeros((N))
    yi = np.zeros((N))
    Ei=Epi=0
    for i in range(N):
        xi[i] = xi[i-1] + np.sin(phi[t_i,i])*l
        yi[i] = yi[i-1] + np.cos(phi[t_i,i])*l
        Ej=[0,0]
        for j in range(i+1): #2. osnovna def.
            Ej[0] += phi_dot[t_i,j]*np.cos(phi[t_i,j])
            Ej[1] += phi_dot[t_i,j]*np.sin(phi[t_i,j])
            Epi += m[i]*g*(L-l*np.cos(phi[t_i,j]))
        Ei += 0.5*m[i]*l**2*(Ej[0]**2 + Ej[1]**2)
    x.append(xi)
    y.append(-yi)
    Ekt.append(Ei) #translacija
    Ep.append(Epi)
E2 = np.array(Ekt)+np.array(Ep)
x2 = np.array(x)
y2 = np.array(y)
Ekt2=np.array(Ekt)/E2[0]
Ep2=np.array(Ep)/E2[0]
E2=E2/E2[0] #normiranje
print('t_post = ', timer()-t2)

##################################
# plotting
def animate(i):
    thisx = np.concatenate(([0], x1[i]))
    thisy = np.concatenate(([0], y1[i]))
    if  i > hist: h = i-hist
    else: h = 0
    line.set_data(thisx, thisy)
    trace.set_data(x1[h:i,-1], y1[h:i,-1])
    time_text.set_text(time_template % (i*dt))
    x_loc.set_data(x1[:i,-1], t[:i])
    y_loc.set_data(t[:i], y1[:i,-1])
    E_t.set_data(t[:i], E1[:i])
    # Ekr_t.set_data(t[:i], Ekr1[:i])
    Ekt_t.set_data(t[:i], Ekt1[:i])
    Ep_t.set_data(t[:i], Ep1[:i])
    ax[0,1].relim()
    ax[0,1].autoscale_view()
    ax[1,0].relim()
    ax[1,0].autoscale_view()
    ax[1,1].relim()
    ax[1,1].autoscale_view()
    return line,time_text,trace,x_loc,y_loc,E_t,Ekt_t,Ep_t
def animate2(i):
    thisx1 = np.concatenate(([0], x1[i]))
    thisy1 = np.concatenate(([0], y1[i]))
    thisx2 = np.concatenate(([0], x2[i]))
    thisy2 = np.concatenate(([0], y2[i]))
    if  i > hist: h = i-hist
    else: h = 0
    history_x1 = x1[h:i,-1]
    history_y1 = y1[h:i,-1]
    history_x2 = x2[h:i,-1]
    history_y2 = y2[h:i,-1]
    line1.set_data(thisx1, thisy1)
    trace1.set_data(history_x1, history_y1)
    line2.set_data(thisx2, thisy2)
    trace2.set_data(history_x2, history_y2)
    time_text.set_text(time_template % (i*dt))
    return line1,line2, trace1,trace2, time_text

fig, ax = plt.subplots(2,2, figsize=(12.8,7.2))
ax[0,0].set_xlim((-L, L))
ax[0,0].set_ylim((-L-l, L))
ax[0,0].set_aspect('equal')
ax[0,0].axis('off')
time_template = 't = %.1f s'
time_text = ax[0,0].text(0, 0, '', transform=ax[0,0].transAxes)
ax[0,0].set_title(f'$N=${N}, $\\varphi_0=${np.rad2deg(phi_01[0]):.0f}$^\circ$, $\\zeta=${d}')
ax[0,1].set_title('$y(t)$')
ax[1,0].set_title('$x(t)$')
ax[1,1].set_title('$E(t)$')
ax[0,1].set_xlabel('$t$')
ax[0,1].set_ylabel('$y$')
ax[1,0].set_xlabel('$x$')
ax[1,0].set_ylabel('$t$')
ax[1,1].set_xlabel('$t$')
ax[1,1].set_ylabel('$E/E_0$')
ax[0,1].set_ylim(-L-l,L)
ax[1,0].set_xlim(-L-l,L+l)
ax[1,1].set_ylim(0,max(E1)*1.2)
ax[0,1].set_xticks([])
ax[1,0].set_yticks([])
ax[1,1].set_xticks([])

line, = ax[0,0].plot([], [], 'k.-', lw=2)
trace, = ax[0,0].plot([], [], 'k-', lw=1, alpha=0.3)
y_loc, = ax[0,1].plot([], [], 'k-',lw=1)
x_loc, = ax[1,0].plot([], [], 'k-',lw=1)
E_t, = ax[1,1].plot([],[], 'r-',label='$E$')
# Ekr_t, = ax[1,1].plot([],[], 'C1-',label='$E_{r}$')
Ekt_t, = ax[1,1].plot([],[], 'b-',label='$E_{t}$')
Ep_t, = ax[1,1].plot([],[], 'g-',label='$E_p$')
ax[1,1].legend(loc='upper right',ncol=4)

ani = animation.FuncAnimation(
    fig, animate, len(t), interval=dt, blit=True)
plt.tight_layout()
plt.show()

# dva plota
fig = plt.figure(figsize=(12.8,7.2)) #HD
ax = fig.add_subplot(autoscale_on=False, xlim=(-L-l, L+l), ylim=(-L-l, L), aspect='equal')
ax.set_aspect('equal')
ax.set_xlabel('$x$ [m]')
ax.set_ylabel('$y$ [m]')
time_template = 't = %.1f s'
time_text = ax.text(0.05, 0.05, '', transform=ax.transAxes)

ax.set_title(f'$N=${N}, $\\zeta={d}$')
line1, = ax.plot([], [], 'b.-', lw=2, label=f'{np.rad2deg(phi_01[0]):.0f}$^\circ$')
trace1, = ax.plot([], [], 'b-', lw=1, alpha=0.5)
line2, = ax.plot([], [], 'g.-', lw=2, label=f'{np.rad2deg(phi_02[0]):.0f}$^\circ$')
trace2, = ax.plot([], [], 'g-', lw=1, alpha=0.5)
ani = animation.FuncAnimation(
    fig, animate2, len(t), interval=dt, blit=True)
ax.legend(title='$\\varphi_0$',loc='upper right')

plt.tight_layout()
#ani.save('./pendulum.gif',writer=animation.PillowWriter(fps=60))#počasi, velik file
#ani.save(f'./pendulum_N{N}_phi{np.rad2deg(phi_0[0]):.0f}.mp4',writer=animation.FFMpegWriter(fps=60))
plt.show()

#energije
plt.plot(t, E1, 'r-',label='$E$')
plt.plot(t, E2, 'r--')
# plt.plot(t, Ekr1, 'b-',label='$E_r$')
# plt.plot(t, Ekr2, 'b--')
plt.plot(t, Ekt1, 'C1-',label='$E_t$')
plt.plot(t, Ekt2, 'C1--')
plt.plot(t, Ep1, 'g-',label='$E_p$')
plt.plot(t, Ep2, 'g--')
legend = plt.legend(handles=[plt.Line2D([], [], color='k', linestyle='-', label=f'{np.rad2deg(phi_01[0]):.0f}$^\circ$'),
                             plt.Line2D([], [], color='k', linestyle='--', label=f'{np.rad2deg(phi_02[0]):.0f}$^\circ$')],
                     loc='upper center', title='$\\varphi_0$',ncol=2)
plt.gca().add_artist(legend)
plt.legend(loc='upper right')
plt.xlabel('$t$ [s]')
plt.ylabel('$E/E_0$')
plt.title(f'$N=${N}, $\\Delta t=${dt}, $\zeta=${d}')
plt.tight_layout()
plt.show()