# -*- coding: utf-8 -*-
# gibanje neraztegljive vrvice
# algoritem po predavanjih Modelska nalaiza II, Neraztegljiva vrvica, 11.pred. 16.5.2024

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import scipy.fft as fft
from timeit import default_timer as timer
from matplotlib.animation import FuncAnimation

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

Ns = 10 #diskretizacija vrvi
ds = 1/Ns
s=np.arange(0,1,ds)

T = 1 #čas simualcije, s
Nt = 5
dt = T/Nt
stabilnost=dt/ds
print('stabilnost:',stabilnost,dt,ds)

#ZP
phi_0 = np.ones((Ns)) * np.deg2rad(45) #zač konst naklon, rad

#inicializacija
phi = phi_0 # phi^n
phi_old = np.copy(phi) # phi^n-1
phi_new = np.zeros_like(phi) # phi^n+1
F = np.zeros_like(phi)

#zanka po času
t1=timer()
sile = []
zasuki = [phi_0]
for t in range(Nt):
    #sila
    diag = -(2*np.ones_like(phi[2:]) + 1/4*(phi[2:]-phi[:-2])**2)
    sp=zg=np.ones((diag.size-1))
    desna = (ds/dt)**2 * (phi_old[1:-1]-phi[1:-1])**2 # = 0 ?!
    F[1:-1] = TDMAsolver(sp,diag,zg,desna)

    #RP sila
    F[0] = F[1] + ds*np.sin(phi[0])
    F[-1] = 0

    #zasuk phi
    phi_new[1:-1] = 2*phi[1:-1] - phi_old[1:-1] + (dt/ds)**2 * (1/2*(F[2:]-F[:-2])*(phi[2:]-phi[:-2]) + F[1:-1]*(phi[2:]-2*phi[1:-1]+phi[:-2]))

    #RP phi
    phi_new[0] = phi_new[2] + 2*ds/(F[0])*np.cos(phi[1])
    phi_new[-1] = 2*phi_new[-2] - phi_new[-3] #rp 2.reda
    #phi_new[-1] = phi_new[-2] #rp 1.red

    zasuki.append(phi_new)
    sile.append(F)

    phi_old = phi
    phi = phi_new
    phi_new = np.zeros_like(phi)
t_sim = timer()-t1
zasuki = np.array(zasuki)
sile = np.array(sile)
print(f'simulation was running: {t_sim} s')

#phi -> x,y
t2=timer()
x = []
y = []
E = [] #energija
for t in range(Nt):
    xt = np.cos(zasuki[t])*ds + s
    yt = np.sin(zasuki[t])*ds + s
    Et = 0 #TODO
    x.append(xt)
    y.append(-yt)#obrnemo os
# for t in range(Nt):
#     xt = np.zeros((Ns))
#     xt[0] = np.sin(zasuki[t,0])*ds
#     yt = np.zeros((Ns))
#     yt[0] = np.cos(zasuki[t,0])*ds
#     for i in range(Ns):
#         xt[i] = xt[i-1] + np.cos(zasuki[t,i])*ds
#         yt[i] = yt[i-1] + np.sin(zasuki[t,i])*ds
#     Et = 0 #TODO
#     x.append(xt)
#     y.append(-yt)#obrnemo os
t_post=timer()-t2
x=np.array(x)
y=np.array(y)
E=np.array(E)
print(f'post process was running: {t_post} s')

#for t in range(Nt):
    #plt.plot(s,zasuki[t],'.-')
# plt.plot(x[-1],y[-1],'.-')
# plt.title(f'{stabilnost:.5f}')
# plt.show()

def animate(i):
    thisx = np.concatenate(([0], x[i]))
    thisy = np.concatenate(([0], y[i]))

    line.set_data(thisx, thisy)
    time_text.set_text(time_template % (i*dt))
    return (line,time_text)

fig = plt.figure(figsize=(5, 4))
ax = fig.add_subplot(autoscale_on=False, xlim=(-1, 1), ylim=(-1-ds, 1))
ax.set_aspect('equal')
ax.grid()
ax.set_xlabel('$x$ [m]')
ax.set_ylabel('$y$ [m]')
ax.set_title(f'$N=${Ns}')

line, = ax.plot([], [], 'k.-', lw=1)
time_template = 't = %.1f s'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

ani = FuncAnimation(
    fig, animate, Nt, interval=500, blit=True)
plt.tight_layout()
plt.show()