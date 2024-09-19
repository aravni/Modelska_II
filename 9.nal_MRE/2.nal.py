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
    return v_0 * (1+b)*y / (y**2+b**4*x**2)**0.5

def NACA(t,N):
    x = np.linspace(1, 0, N//2, endpoint=True)
    yt = t/50 *(1.457122 * x**0.5 -
                0.624424 * x -
                1.727016 * x**2 +
                1.384087 * x**3 -
                0.489769 * x**4)
    #zlepimo obe polovici
    x = np.concatenate([x, x[-2::-1]])
    y = np.concatenate([yt, -yt[-2::-1]])
    return x,y

def žukovski(phi,A,B,R=1):
    z = A+B*1j+R*np.exp(phi*1j)
    r =  np.array(1/2*(z+1/z))
    x=r.real
    y=r.imag
    x[-1]=x[0]
    y[-1]=y[0]
    return x[::-1],y[::-1]

def žukovski2(phi,m,delta,R):
    theta = R*(m*np.exp(1j*delta) + np.exp(phi*1j))
    z = theta + 1/theta
    x=z.real
    y=z.imag
    x[-1]=x[0]
    y[-1]=y[0]
    return x[::-1],y[::-1]

def profil(x_nodes, y_nodes):
    '''tvori območje z paneli na vozliščnih koordinatah'''

    # get middle points from nodes
    N=len(x_nodes)
    x_p = []
    y_p = []
    l_p = []
    phi_p = []
    for i in range(N-1):
        # xi = (x_nodes[i-1]+x_nodes[i])/2
        # yi = (y_nodes[i-1]+y_nodes[i])/2
        # li = ((x_nodes[i-1]-x_nodes[i])**2 + (y_nodes[i-1]-y_nodes[i])**2)**0.5
        # #phi_i = np.arctan2(y_nodes[i-1]-y_nodes[i], x_nodes[i-1]-x_nodes[i]) #normale ven
        # phi_i = np.arctan2(y_nodes[i]-y_nodes[i-1], x_nodes[i]-x_nodes[i-1]) #?!
        xi = (x_nodes[i]+x_nodes[i+1])/2
        yi = (y_nodes[i]+y_nodes[i+1])/2
        li = ((x_nodes[i+1]-x_nodes[i])**2 + (y_nodes[i+1]-y_nodes[i])**2)**0.5
        phi_i = np.arctan2(y_nodes[i+1]-y_nodes[i], x_nodes[i+1]-x_nodes[i]) #štetje nasproti ure
        #phi_i = np.arctan2(y_nodes[-i-1]-y_nodes[-i], x_nodes[-i-1]-x_nodes[-i]) #normale ven! v smeri ure?
        #print(i,xi,yi,li,phi_i)
        # if yi == 0.: #singularno
        #     yi = 1e-10
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
    b_eq = np.zeros((N))
    for i in range(N):
        b_eq[i] = -v_0[0] * np.sin(phi_p[i]) + v_0[1] * np.cos(phi_p[i])
        #print(b[i])
        for j in range(N):
            # xi, yi = R_gl(phi_p[j]) @ np.array([x_p[i]-x_p[j],y_p[i]-y_p[j]]) #lega panela i v lokalnem KS panela j
            # #print(i,j,':',xi,yi,phi_p[j])
            # if yi == 0.:
            #     yi = 1e-10
            # v_lok_j = np.array([v_parallel(xi, yi, l=l_p[j]), v_perpendicular(xi, yi, l=l_p[j])]) #hitrosti na mestu i-panela v lokalnem KS j-panela
            # # v_lok_j = np.array([0,1])
            # #print(xi,yi,':',v_lok_j)
            # v_gl = R_lok(phi_p[j]) @ v_lok_j #hitrosti v globalnem KS
            # #print(v_gl)
            # v_lok_i = R_gl(phi_p[i]) @ v_gl #hitrosti v lokalnem KS i-panela
            # #print(v_lok_i)

            #github
            # rotate coordinates into CM frame
            xp = ((x_p[j] - x_p[i]) * np.cos(phi_p[i]) +
                (y_p[j] - y_p[i]) * np.sin(phi_p[i]))
            yp = (-(x_p[j] - x_p[i]) * np.sin(phi_p[i]) +
                (y_p[j] - y_p[i]) * np.cos(phi_p[i]))
            # if yp == 0.:
            #     yp = 1e-10
            # calculate vp and vt from formula
            v = np.array([v_perpendicular(xp, yp, l_p[i]),  v_parallel(xp, yp, l_p[i])])

            # rotate vp/vt into vp in j frame
            A[i,j] += (v[0] * np.cos(phi_p[i] - phi_p[j]) +
                v[1] * np.sin(phi_p[i] - phi_p[j]))
            #print(A[i,j], x_p[j], y_p[j], phi_p[j], l_p[j], v, xp,yp)
        #print(x_p[i], y_p[i], phi_p[i], l_p[i])
    sigma = sc.linalg.solve(A, b_eq, overwrite_a=True, overwrite_b=True)
    c = 1/2 * np.sum(sigma*l_p)
    #print('sigma',sigma)

    M=20 #to ja drago!
    meja=3
    x = np.linspace(-meja,meja,M)
    y = np.linspace(-meja,meja,M)
    X, Y = np.meshgrid(x, y)

    u = np.zeros((M,M))
    v = np.zeros((M,M))
    for i in range(M):
        for j in range(M):
            sx = v_0[0]
            sy = v_0[1]
            for n in range(N):
                # xi, yi = R_gl(phi_p[n]) @ np.array([x[i]-x_p[n],y[j]-y_p[n]]) #lega točke (i,j) v lokalnem KS n-panela
                # if yi == 0.:
                #     yi = 1e-10
                # v_lok = np.array([v_parallel(xi, yi, l=l_p[n]), v_perpendicular(xi, yi, l=l_p[n])]) # hitrosti v lokalnem KS n-panela
                # v_gl = R_lok(phi_p[n]) @ v_lok #hitrosti v globalnem KS
                # #print(x[i],y[j],xi,yi,':', v_lok, v_gl)
                # u[j,i] += v_gl[0] * sigma[n]
                # v[j,i] += v_gl[1] * sigma[n]

                #github
                # rotate coordinates into CM frame
                # xp = ((x[i] - x_p[n]) * np.cos(phi_p[n]) +
                #     (y[j] - y_p[n]) * np.sin(phi_p[n]))
                # yp = (-(x[i] - x_p[n]) * np.sin(phi_p[n]) +
                #     (y[j] - y_p[n]) * np.cos(phi_p[n]))
                xp = ((X[i,j] - x_p[n]) * np.cos(phi_p[n]) +
                    (Y[i,j] - y_p[n]) * np.sin(phi_p[n]))
                yp = (-(X[i,j] - x_p[n]) * np.sin(phi_p[n]) +
                    (Y[i,j] - y_p[n]) * np.cos(phi_p[n])) #zakaj je pa to OK?! -> pravilno plotanje !!!
                # contibution in CM frame from formula
                v_ = sigma[n] * np.array([v_perpendicular(xp, yp, l_p[n]),  v_parallel(xp, yp, l_p[n])])
                # rotate vp/vt into LAB frame
                vy = v_[0] * np.cos(phi_p[n]) + v_[1] * np.sin(phi_p[n])
                vx = v_[1] * np.cos(phi_p[n]) - v_[0] * np.sin(phi_p[n])
                #print(v_,vx,vy, u[i,j])
                # add contributions
                sx += vx
                sy += vy
            # if (x[i]/np.cos(np.arctan(y[j]/x[i])))**2 + (y[j]/(b*np.sin(np.arctan(y[j]/x[i]))))**2 > 1: #le za elipse!!
            #     u[i,j] = sx
            #     v[i,j] = sy
            u[i,j] = sx
            v[i,j] = sy
            #print(x[i],y[j],':', u[i,j],v[i,j])

    fig,ax =plt.subplots(1, figsize=(8,7))
    ax.plot(np.append(x_nodes, x_nodes[0]), np.append(y_nodes, y_nodes[0]), '.-',c='k')

    # #vec = ax.quiver(x,y, u.T, v.T, (u**2+v**2)**0.5) #večji scale, manjši vektor!
    # vec = ax.quiver(X,Y, u, v, (u**2+v**2)**0.5) #meshgrid!
    # fig.colorbar(vec, ax=ax, label='|$\mathbf{v}$|')

    strm = ax.streamplot(X,Y, u,v, broken_streamlines=False, color=(u**2+v**2)**0.5, cmap='turbo', density=1)
    fig.colorbar(strm.lines, label='|$\mathbf{v}$|')
    for art in ax.get_children(): #remove arrows
        if not isinstance(art, matplotlib.patches.FancyArrowPatch):
            continue
        art.remove()        # Method 1
        #art.set_alpha(0)  # Method 2

    #plt.title(f'$N=${N}')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.axis('equal')

#hitrost v neskončnosti
alpha = np.deg2rad(0) #angle of attack
v_0 = np.array([np.cos(alpha), np.sin(alpha)]) *14 #pod kotom
print(v_0)
N = 50

#elipsa
b = 1.5
phi = np.linspace(0,2*np.pi,N+1,endpoint=True)
x,y = elipsoid(phi, b)
title = f'$N=${N}, $b=${b}, $\\alpha=${np.rad2deg(alpha):n}˚'

#NACA
# t=15
# x,y = NACA(t=t, N=N)
# title = f'$N=${N}, $t=${t}'

#Žukovski
# A,B = -0.2, 0.3 #airfoil
# R = ((1+abs(A))**2+(A**2+B**2))**0.5
# #A,B,R = -0.2, 0.1, 1.25 #airfoil
# A,B,R = 0,0,0.5 #elipsa
# phi = np.linspace(0,2*np.pi,N+1,endpoint=True)
# x,y = žukovski(phi, A,B,R) #A<0 - nos levo, B>0 - krivljenje, R=1 - ravna črta, navodila R<1
# title = f'$N=${N}, $A=${A}, $B=${B}, $R=${R:.3f}'
# plt.plot(x,y,'.-')
# plt.axis('equal')
# plt.show()

# m=0.51
# delta = np.pi*0.66 #pi=brez krivljena, delta>pi/2 nos je levo
# R = ((1+m)**2+(m*np.tan(delta))**2)**0.5
# #R=1.510112886
# print(R)
# x,y = žukovski2(phi, m, delta, R=R) #polarni popis
# plt.plot(x,y,'.-')
# plt.grid()
# plt.title(title)
# plt.show()

#clark-Y
# title = 'Clark-Y airfoil'
# x1,y1 = np.loadtxt('clark_y_up.txt', unpack=True)
# x2,y2 = np.loadtxt('clark_y_low.txt', unpack=True)
# n=3 #korak točk
# x2 = np.delete(x2, 0)
# y2 = np.delete(y2, 0)
# x = np.append(x1[::-n],x2[n-1::n])
# y = np.append(y1[::-n],y2[n-1::n])
# y[-1] = y[0]
# x=x[::-1]
# y=y[::-1]
# # plt.plot(x,y,'.-')
# plt.show()

#EPPLER 435
title='EPPLER 435 MRE'
x1,y1 = np.loadtxt('eppler435_up.dat', unpack=True)
x2,y2 = np.loadtxt('eppler435_low.dat', unpack=True)
n=5 #korak točk
x = np.append(x1[::-n],x2[n+1::n])
y = np.append(y1[::-n],y2[n+1::n])

######### solve
profil(x,y)
plt.title(title)
plt.axis('equal')
plt.tight_layout()
#plt.savefig(f'./b{b}_N{N}_alpha{np.rad2deg(alpha):n}.pdf')#elipsa
#plt.savefig(f'./NACA_vec_t{t}_N{N}_alpha{np.rad2deg(alpha):n}.pdf')#NACA
#plt.savefig(f'./zukovski_A{A}_B{B}_R{R:.3f}_N{N}_alpha{np.rad2deg(alpha):n}.pdf')#žukovski
#plt.savefig(f'./clarkY_alpha{np.rad2deg(alpha):n}.pdf')
plt.show()