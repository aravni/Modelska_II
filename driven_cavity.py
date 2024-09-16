# coding: utf-8
# reševanje problema gnane kotanje s SIMPLE in PISO algoritmom

#import SIMPLE
from SIMPLE import *
#from PISO import PISO_old, PISO
from timeit import default_timer as timer
import numpy as np

N = 10 #št. diskretizacijskih točk
FD_kind = 'hybrid' #diskretizacijska shema ['central','upwind', 'hybrid']
Re = 10 #Reynolds
alpha_p = 1 #relaksacija tlaka
alpha_vel = 1 #relaksacija hitrosti
N_max_iters = 100 #max št iteracij

tol = 1e-3 #toleranca ohranitve mase
v_lid = 1 #hitrost na zg robu -> vpliva na Re -> če = 1 je brezdim!
b = 1 #debelina KV
p0 = 1 #začetni približek tlaka
RP = (1,1,0,0)

#PISO(N, FD_kind, Re, alpha_p, alpha_vel, N_max_iters, v_lid, rho, tol)
#F, time, continuity, iter = PISO(N, FD_kind, Re, alpha_p, alpha_vel, N_max_iters, tol, b, piso=True, p_0=p0, RP=RP)
rho=1
F, time, continuity, iter = SIMPLE(N, FD_kind, Re, alpha_p, alpha_vel, N_max_iters, v_lid, rho, tol, b)

w = input('Write data to file? y/n')

if w.upper() == 'Y':
    with open('./results_SIMPLE.txt','a') as d:
        d.write(f'\n{Re}\t{FD_kind}\t{N*N}\t{iter}\t{time}\t{F[0]}\t{F[1]}\t{F[2]}\t{F[3]}\t{continuity[-1]:E}\t{alpha_p}\t{alpha_vel}')
    print('Data written to file!')
else: print('Data not written!')