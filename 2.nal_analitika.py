# -*- coding: utf-8 -*-
# obtekanje teles

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib
import scipy as sc
from matplotlib.cm import ScalarMappable

rc({"text.usetex": True,
    'family': 'serif', 'serif': ['Computer Modern']})
plt.rcParams.update({'font.size': 10})
def pause():
    input("Press the <ENTER> key to continue...")

def R_gl(phi):
    '''rotacija: globalni -> lokalni'''
    return np.array([[np.cos(phi), np.sin(phi)],[-np.sin(phi), np.cos(phi)]])
def R_lok(phi):
    '''rotacija: lokalni -> globalni'''
    return np.array([[np.cos(phi), -np.sin(phi)],[np.sin(phi), np.cos(phi)]])

def elipsoid(phi,b):
    return np.cos(phi), b*np.sin(phi)

def v_parallel(x,y,l): #y=0 singularnost!!!
    if y == 0 and (x-l/2 > 0 or x+l/2 < 0):
        return 1 / (2 * np.pi) * np.log(np.abs((x+l/2) / (x-l/2)))
    elif y == 0 and x-l/2 < 0 and x+l/2 > 0:
        return 1 / (2 * np.pi) * np.log(np.abs((x+l/2) / (x-l/2)))
    elif x-l/2 == 0:
        return 1 / (4 * np.pi) * ((x+l/2)/y)**2
    elif x+l/2 == 0:
        return -1 / (4 * np.pi) * ((x-l/2)/y)**2
    else:
        return 1/(4*np.pi)* np.log(((x+l/2)**2+y**2) / ((x-l/2)**2+y**2))
def v_perpendicular(x,y,l):
    if y == 0 and (x-l/2 > 0 or x+l/2 < 0):
        return 0
    elif y == 0 and x-l/2 < 0 and x+l/2 > 0:
        return 1/2
    elif x-l/2 == 0:
        return 1 / (2 * np.pi) * np.arctan(x+l/2 / y)
    elif x+l/2 == 0:
        return 1 / (2 * np.pi) * np.arctan(x-l/2 / y)
    else:
        return 1/(2*np.pi)* (np.arctan((x+l/2)/y) - np.arctan((x-l/2)/y))
    
def v_parallel_analitical(x,y,b,v_0):
    '''analitičen tok ob površini elipse'''
    l=np.zeros_like(x)
    for i in range(len(x)):
        dx = np.diff(x[:i+1])
        dy = np.diff(y[:i+1])
        segment_lengths = np.sqrt(dx**2 + dy**2)
        l[i] = np.sum(segment_lengths)
    v = v_0 * (1+b)*y / (y**2+b**4*x**2)**0.5
    return l, v

def profil(x_nodes, y_nodes, b_elipse):
    '''tvori območje z paneli na vozliščnih koordinatah'''

    # get middle points from nodes
    N=len(x_nodes)
    x_p = []
    y_p = []
    l_p = []
    phi_p = []
    for i in range(N-1):
        xi = (x_nodes[i]+x_nodes[i+1])/2
        yi = (y_nodes[i]+y_nodes[i+1])/2
        li = ((x_nodes[i+1]-x_nodes[i])**2 + (y_nodes[i+1]-y_nodes[i])**2)**0.5
        phi_i = np.arctan2(y_nodes[i+1]-y_nodes[i], x_nodes[i+1]-x_nodes[i]) #štetje nasproti ure
        x_p.append(xi)
        y_p.append(yi)
        l_p.append(li)
        phi_p.append(phi_i)
    x_p = np.array(x_p)
    y_p = np.array(y_p)
    l_p = np.array(l_p)
    phi_p = np.array(phi_p)

    N=len(x_p) #št. panelov
    A = np.zeros((N,N))
    b = np.zeros((N))
    for i in range(N):
        b[i] = -v_0[0] * np.sin(phi_p[i]) + v_0[1] * np.cos(phi_p[i])
        for j in range(N):
            # rotate coordinates into CM frame
            xp = ((x_p[j] - x_p[i]) * np.cos(phi_p[i]) +
                (y_p[j] - y_p[i]) * np.sin(phi_p[i]))
            yp = (-(x_p[j] - x_p[i]) * np.sin(phi_p[i]) +
                (y_p[j] - y_p[i]) * np.cos(phi_p[i]))
            # calculate vp and vt from formula
            v = np.array([v_perpendicular(xp, yp, l_p[i]),  v_parallel(xp, yp, l_p[i])])
            # rotate vp/vt into vp in j frame
            A[i,j] += (v[0] * np.cos(phi_p[i] - phi_p[j]) +
                v[1] * np.sin(phi_p[i] - phi_p[j]))
    sigma = sc.linalg.solve(A, b, overwrite_a=True, overwrite_b=True)
    c = 1/2 * np.sum(sigma*l_p)

    M=N #računske točke na robu (panelih)
    phi = np.linspace(0,2*np.pi,M+1,endpoint=True)
    x,y = elipsoid(phi, b=b_elipse)

    vt = np.zeros((M))
    l = np.zeros((M))
    for i in range(M):
        st = 0#v_0[0]
        for n in range(N):
            xp = ((x[i] - x_p[n]) * np.cos(phi_p[n]) +
                (y[i] - y_p[n]) * np.sin(phi_p[n])) *(1+1e-10)
            yp = (-(x[i] - x_p[n]) * np.sin(phi_p[n]) +
                (y[i] - y_p[n]) * np.cos(phi_p[n])) *(1+1e-10)
            # contibution in CM frame from formula
            v_ = sigma[n] * np.array([v_perpendicular(xp, yp, l_p[n]),  v_parallel(xp, yp, l_p[n])])
            # rotate vp/vt into LAB frame
            vy = v_[0] * np.cos(phi_p[n]) + v_[1] * np.sin(phi_p[n])
            vx = v_[1] * np.cos(phi_p[n]) - v_[0] * np.sin(phi_p[n])
            # rotate into i frame - vt parallel
            vt_ = vx * np.cos(phi_p[n]) + vy * np.sin(phi_p[n])
            # add contributions
            st += -vt_
        vt[i] = st
        dx = np.diff(x[:i+1])
        dy = np.diff(y[:i+1])
        segment_lengths = np.sqrt(dx**2 + dy**2)
        l[i] = np.sum(segment_lengths)
        #print(i,':',s_parallel)
    return l, vt

#hitrost v neskončnosti
alpha = np.deg2rad(0) #angle of attack
v_0 = [np.cos(alpha), np.sin(alpha)] #pod kotom

#analitika
b=1
phi = np.linspace(0,2*np.pi,500+1,endpoint=True)
x,y = elipsoid(phi, b=b)
l_analit, v_analit = v_parallel_analitical(x,y,b=b,v_0=v_0[0])
plt.plot(l_analit, v_analit,'k-',lw=5,label='analitično')

res=[]
N_list = np.arange(10,101,1,dtype=int)#residuali
N_list = np.arange(20,101,20,dtype=int)#plot
for N in N_list:
    phi = np.linspace(0,2*np.pi,N+1,endpoint=True)
    x,y = elipsoid(phi, b=b)
    l, vt = profil(x,y,b)
    l_a, v_a = v_parallel_analitical(x,y,b=b,v_0=v_0[0])
    res.append(np.sum(np.abs(vt-v_a[:len(vt)])))
    plt.plot(l, vt, '.-', label=f'{N}')

plt.grid()
plt.legend(title='št.panelov')
plt.title(f'tangencialna hitrost ob površini elipse $b=${b}')
plt.xlabel('$l$')
plt.ylabel('$v_{\parallel}$')
#plt.savefig(f'./vt_b{b}.pdf')
plt.show()

plt.plot(N_list, res,'.-')
plt.grid()
plt.xlabel('$N$')
plt.ylabel('$|v_{parallel}-v_{analitično}|$')
plt.title(f'residual glede na število panelov, $b=${b}')
#plt.savefig(f'./residual_b{b}.pdf')
plt.show()
print(min(res))