#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Zaključna naloga: upogib nosilca, statika
# Euler-Bernoulli nosilec, PLASTIČNOST
# viri:
# [1] E. A. de Souza Neto, D. Peric, D. R. J. Owen - Computational Methods for Plasticity-Theory and Applications - Wiley (2008)
# [2] https://vulcanhammer.net/wp-content/uploads/2017/01/finite-elements-in-plasticity-theory-and-practice.pdf

import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from scipy.sparse.linalg import spsolve
from timeit import default_timer as timer
from alignAxes import *
from colored_line import colored_line
import findiff

plt.style.use('aljaz')

class solver:
    def __init__(self, Ne, Nh, N_iter, tol):
        self.Ne=Ne #št elementov
        self.Nv=Ne+1 #št vozlišč
        self.Nh=Nh #št točk po višini prereza
        self.N_iter = N_iter
        self.tol=tol

        self.load_incr = 1 #faktor obremenitve
        self.Nps = 2 #št stopenj prostosti na vozlišče
        self.KE = 'EB' #izbira KE (EB=nosilec Euler-Bernoulli, TG1,TG2,TM-oblike za Timošenko)

        # snovne lastnosti
        self.E, self.R_p, self.H0 = 1, 1, 0.5 #test
        # self.R_p = 300 # meja plastičnosti, MPa
        # ρ = 1 # kg/mm^3
        # self.E = 189126 # MPa
        self.ν = 0.3 # Poisson
        # Et = 0.003*E # konst. tangentni modul v plastičnem, MPa (np.nan za nelinearno)
        self.H=self.H0 #= Et/(1-Et/E) # linearno utrjevanje
        self.G = self.E/(2*(1+self.ν)) # strižni modul v elast. območju, MPa
        self.K=self.E/(3*(1-2*self.ν)) # bulk modul, MPa
    
        # obremenitve
        self.RP = [2,0] # vpetje (0-prosto, 1-fiksno, 2-konzolno) [levo, desno]
        self.RP_F = [0,0,1] # točkovna sila F0 (0-neobremenjeno, 1-obremenjeno) [levo, sredina, desno]
        self.RP_M = [0,0,0] # upogibni moment M0 (0-neobremenjeno, 1-obremenjeno) [levo, sredina, desno]
        self.F0 = 0.1 # točkovna sila - navzdol, N
        self.M0 = 500*1e3 # upogibni moment - smer ure, N mm
        self.q0 = 0 # linijska obremenitev, N/mm
    
        # geometrijske lastnosti
        self.L = 2 # dolžina nosilca, mm
        # 1. pravokotnik: a-širina, h-višina
        self.a, self.h = 1, 1 #mm
        κ = 10*(1+self.ν)/(12+11*self.ν) #5/6 #koef. oblike prereza
        self.A = self.a*self.h #mm^2
        self.I = 1/12 * self.a*self.h**3 #mm^4
        # 2. krog: r-polmer
        # r = 100 #mm
        # κ = 6*(1+self.ν)/(7+6*self.ν) #9/10 #koef. oblike prereza
        # A = np.pi*r**2 #mm^2
        # I = 1/64 * np.pi*r**4 #mm^4
    
        self.hx = self.L/(self.Nv-1) # dolžina KE
        self.hy = self.h/(self.Nh-1) # višina deljenja prereza
        As = κ*self.A # vpliv oblike prereza na velikost max strižnih nap., upoštevajoč konst. porazd. stižnih nap.
        M_pl = self.R_p * self.a * self.h**2/4 # notranji moment pri meji plastičnosti
        A_layer = self.a * self.hy # površina sloja pri višini z

        # TODO preveri ali so dimenzije ok
        self.ε = np.zeros((Ne,Nh,3,3))
        self.ε_e = np.zeros((Ne,Nh,3,3))
        self.ε_p_eq = np.zeros((Ne,Nh))

        self.Φ = np.zeros((Ne,Nh))
        self.σ_y = np.zeros((Ne,Nh))

        self.p = np.zeros((Ne,Nh))
        self.s = np.zeros((Ne,Nh,3,3))
        self.q = np.zeros((Ne,Nh))
        self.σ = np.zeros((Ne,Nh,3,3))

        self.f_ext = np.zeros((self.Nv*self.Nps))
        self.f_int = np.zeros((self.Nv*self.Nps))
        self.r = np.zeros((self.Nv*self.Nps))
        self.M = np.zeros(Ne)
        self.T = np.zeros(Ne)

        self.u = np.zeros((self.Nv*self.Nps))
        self.U = np.zeros((self.Nv))
        self.Θ = np.zeros((self.Nv))
        self.u_integ = np.zeros((Ne))

        self.identity = np.eye(3)
        if Nh == 1: self.y_list = np.array([self.h/2])
        else:
            # self.y_list = np.linspace(-self.h/2., self.h/2., Nh) #višina prereza
            self.y_list = np.linspace(0, self.h/2., Nh) #višina prereza SYM

    def upogibnica(self, n=None):
        "analitična rešitev konzolno vpetega Euler Bernoulli nosilca"
        if n is None: x = np.linspace(0,self.L,self.Nv) #vozlišča
        else: x = np.linspace(0,self.L,n) #n točk
        w = self.load_incr * self.F0*x**2/(6*self.E*self.I)*(3*self.L-x) #+ M0*x**2/(2*self.E*self.I)
        return w

    def shape_f(self,xi):
        "calculates shape functions and their derivatives for a 2-node beam element of length L at a given local coordinate position xi (-1 <= xi <= 1)"
        L=self.hx #dolžina KE
        N = np.array([
            0.25*(2-3*xi+xi**3),
            L/8*(1-xi-xi**2+xi**3),
            0.25*(2+3*xi-xi**3),
            L/8*(-1-xi+xi**2+xi**3)
        ])
        dNdx = np.array([
            3*(xi**2-1)/(2*L),
            (-1-2*xi+3*xi**2)/4,
            3*(1-xi**2)/(2*L),
            (-1+2*xi+3*xi**2)/4
        ])
        B = np.array([
            -6*xi/L**2,
            -(3*xi-1)/L,
            6*xi/L**2,
            -(3*xi+1)/L
        ]) #d²N/dx²
        return N, dNdx, B
    
    # diskretizacija
    #############################################################################

    def K_e_eb(self, Le):
        """togostna matrika za KE Euler-Bernoulli nosilec"""
        K_lok = self.E*self.I/Le**3*np.array([[12,   6*Le,   -12,   6*Le],
                                    [6*Le, 4*Le**2,-6*Le, 2*Le**2],
                                    [-12, -6*Le,    12,  -6*Le],
                                    [6*Le, 2*Le**2,-6*Le, 4*Le**2]])
        return K_lok
    def K_e_tg2(self, Le):
        """togostna matrika za KE Timošenkov nosilec Gaussova integracija v 2 točkah
        - pretoge rešitve po [https://vulcanhammer.net/wp-content/uploads/2017/01/finite-elements-in-plasticity-theory-and-practice.pdf]"""
        M1 = self.E*self.I/Le*np.array([[0, 0, 0, 0],
                            [0, 1, 0,-1],
                            [0, 0, 0, 0],
                            [0,-1, 0, 1]])
        M2 = self.G*self.As/Le*np.array([[1,     Le/2,   -1,     Le/2],
                            [Le/2,  Le**2/3,-Le/2,  Le**2/6],
                            [-1,   -Le/2,    1,    -Le/2],
                            [Le/2,  Le**2/6,-Le/2,  Le**2/3]])
        return M1+M2
    def K_e_tg1(self, Le):
        """togostna matrika za KE Timošenkov nosilec Gaussova integracija v 1 točki"""
        M1 = self.E*self.I/Le*np.array([[0, 0, 0, 0],
                            [0, 1, 0,-1],
                            [0, 0, 0, 0],
                            [0,-1, 0, 1]])
        M2 = self.G*self.As/Le*np.array([[1,     Le/2,   -1,     Le/2],
                            [Le/2,  Le**2/4,-Le/2,  Le**2/4],
                            [-1,   -Le/2,    1,    -Le/2],
                            [Le/2,  Le**2/4,-Le/2,  Le**2/4]])
        return M1+M2
    def K_e_tm(self, Le):
        """togostna matrika za KE Timošenkov nosilec modificirana oblika
        - kubična polinomska aproksimacija povesa in zasuka -> OK!"""
        C = (12*self.E*self.I)/(self.As*self.G*Le**2)
        K_lok = self.E*self.I/Le**3/(1+C)*np.array([[12,   6*Le,       -12,   6*Le],
                                        [6*Le, (4+C)*Le**2,-6*Le, (2-C)*Le**2],
                                        [-12, -6*Le,        12,  -6*Le],
                                        [6*Le, (2-C)*Le**2,-6*Le, (4+C)*Le**2]])
        return K_lok

    def discretization(self):
        "diskretizacija, sestava matrik, RP"
        # vozlišča in elementi
        self.vozlisca = np.transpose([np.linspace(0,self.L,self.Nv), np.zeros(self.Nv)])
        ind = np.linspace(0,self.Nv-1,self.Nv, dtype=int)
        elementi = np.transpose([ind[:-1], ind[1:]]).astype(int)
    
        # sestav globalne masne in togostne matrike
        M_glob = np.zeros((self.Nv*self.Nps,self.Nv*self.Nps))
        self.K_glob = np.zeros((self.Nv*self.Nps,self.Nv*self.Nps))
        self.F_glob = np.zeros(self.Nv*self.Nps) # desna stran - vozliščne obremenitve
    
        def distance(a,b): return ((a[0]-b[0])**2 + (a[1]-b[1])**2)**0.5
        for element in elementi:
            Le = distance(self.vozlisca[element[0]], self.vozlisca[element[1]]) # dolžina elementa (Le = hx)
            # priševanje prispevkov posameznih elementov
            ind = np.array([[self.Nps*element[0],self.Nps*element[0]+1,self.Nps*element[1],self.Nps*element[1]+1]])
            if self.KE == 'EB': self.K_glob[ind.T,ind] += sol.K_e_eb(Le)
            elif self.KE == 'TG1': self.K_glob[ind.T,ind] += sol.K_e_tg1(Le)
            elif self.KE == 'TG2': self.K_glob[ind.T,ind] += sol.K_e_tg2(Le)
            elif self.KE == 'TM': self.K_glob[ind.T,ind] += sol.K_e_tm(Le)
            else: raise Exception('Napačen tip KE!')
            self.F_glob[ind[0,[0,2]]] += -self.q0*Le/2
    
        # robni pogoji
        if self.RP[0] == 1 or self.RP[0] == 2: #prvo vozlisce je fiksno
            i = 0
            self.K_glob[i,:] = 0
            self.K_glob[i,i] = 1
            self.F_glob[i] = 0
        if self.RP[0] == 2: #prvo vozlisce konzola
            i = 1
            self.K_glob[i,:] = 0
            self.K_glob[i,i] = 1
            self.F_glob[i] = 0
        if self.RP[1] == 1 or self.RP[1] == 2: #zadnje vozlisce je fiksno
            i = -2
            self.K_glob[i,:] = 0
            self.K_glob[i,i] = 1
            self.F_glob[i] = 0
        if self.RP[1] == 2: #zadnje vozlisce konzola
            i = -1
            self.K_glob[i,:] = 0
            self.K_glob[i,i] = 1
            self.F_glob[i] = 0
    
        # obremenitve
        if self.RP_F[0] == 1: #prvo vozlisce sila F0
            self.F_glob[0] += -self.F0
        if self.RP_F[1] == 1: #središčno vozlisce sila F0
            self.F_glob[(self.Nv*self.Nps)//2] += -self.F0
        if self.RP_F[2] == 1: #zadnje vozlisce sila F0
            self.F_glob[-2] += -self.F0
    
        if self.RP_M[0] == 1: #prvo vozlisce moment
            self.F_glob[1] += -self.M0
        if self.RP_M[1] == 1: #središčno vozlisce moment
            self.F_glob[(self.Nv*self.Nps)//2+1] += -self.M0
        if self.RP_M[2] == 1: #zadnje vozlisce moment
            self.F_glob[-1] += -self.M0

        # geometrija
        self.integ_tocke = np.transpose([sol.vozlisca[:-1,0] + self.hx/2, sol.vozlisca[:-1,1]]) #sredina KE
        
        self.X_defor, self.Y_defor = np.empty((self.Nv,Nh)), np.empty((self.Nv,Nh))
        # for i in range(Nh): X[:,i], Y[:,i] = sol.vozlisca[:,0], sol.vozlisca[:,1] + hy*i - h/2 #cel prerez
        for i in range(Nh): self.X_defor[:,i], self.Y_defor[:,i] = sol.vozlisca[:,0], sol.vozlisca[:,1] + self.hy*i #sym

    # iteracijski postopek
    #############################################################################

    def local1d(self, Δε):
        "samo sigma_xx - Euler-Bernoulli nosilec"
        self.ε[:,:,0,0] += Δε[:,:,0,0] #eps_xx skupna
        self.ε[:,:,1,1] -= Δε[:,:,0,0]*self.ν #eps_yy
        self.ε[:,:,2,2] -= Δε[:,:,0,0]*self.ν #eps_zz

        ε_e_trial = self.ε_e[:,:,0,0] + Δε[:,:,0,0] #eps_xx
        q_trial = self.E * np.abs(ε_e_trial) #euler bernoulli! sigma_xx

        self.σ_y[:,:] = self.R_p + self.H0*self.ε_p_eq #set linear hardening
        phi_trial = q_trial - self.σ_y

        plast=0
        for i in range(self.Ne):
            for j in range(self.Nh): #višina prereza - integracija
                if phi_trial[i,j] <= 0: #elastično
                    self.Φ[i,j] = phi_trial[i,j]
                    self.ε_e[i,j,0,0] = ε_e_trial[i,j]
                    self.q[i,j] = q_trial[i,j]
                    self.σ[i,j,0,0] = q_trial[i,j]
                else: #plastično
                    plast+=1
                    self.ε_e[i,j,0,0] = ε_e_trial[i,j] * self.R_p / q_trial[i,j] #radial return mapping eps_xx
                    self.ε_e[i,j,1,1] = 0.5*self.ε[i,j,0,0] - self.ε_e[i,j,0,0]*(self.ν+0.5) #eps_yy
                    self.ε_e[i,j,2,2] = self.ε_e[i,j,1,1] #eps_zz

                    self.ε_p_eq[i,j] = abs(self.ε[i,j,0,0] - self.ε_e[i,j,0,0])
                    self.σ_y[i,j] = self.R_p + self.H0*self.ε_p_eq[i,j]
                    self.σ[i,j,0,0] = self.σ_y[i,j] #sigma_xx
                    self.q[i,j] = self.σ[i,j,0,0]
                    self.Φ[i,j] = 0. #for numerical stability
        print(f'plastičnih točk: {plast} / {self.Ne*self.Nh}')

    def local_loop(self, Δε): #TODO test pravilnosti lokalne zanke!!!
        ε_e_trial = self.ε_e + Δε
        ε_e_v_trial = np.trace(ε_e_trial,axis1=2,axis2=3)/3.
        # print(ε_e_v_trial.shape, ε_e_trial.shape)
        ε_e_d_trial = ε_e_trial - ε_e_v_trial[:,:,None,None] * self.identity
        # print(ε_e_d_trial.shape)
        p_trial = self.K*ε_e_v_trial
        s_trial = 2.*self.G*ε_e_d_trial

        q_trial = np.sqrt(3./2. * np.sum(s_trial * s_trial, axis=(2,3)))

        self.σ_y[:,:] = self.R_p + self.H0*self.ε_p_eq #set linear hardening
        phi_trial = q_trial - self.σ_y

        plast=0
        for i in range(self.Ne):
            for j in range(self.Nh): #višina prereza - integracija
                if phi_trial[i,j] <= 0: #elastično
                    self.Φ[i,j] = phi_trial[i,j]
                    self.ε_e[i,j] = ε_e_trial[i,j]
                    self.p[i,j] = p_trial[i,j]
                    self.s[i,j] = s_trial[i,j]
                    self.q[i,j] = q_trial[i,j]
                    # self.ε_p_eq[i,j] += 0 #no plastic increment
                else: #plastično
                    plast+=1
                    # print('plastic', i,j)
                    Δγ = phi_trial[i,j] / (3.*self.G + self.H0) #linear isotropic hardening
                    self.ε_p_eq[i,j] += Δγ #plastic increment
                    self.σ_y[i,j] = self.R_p + self.H0*self.ε_p_eq[i,j]

                    self.p[i,j] = p_trial[i,j] #no change in volumetric stress
                    self.s[i,j] = (1. - 3.*Δγ*self.G / q_trial[i,j]) * s_trial[i,j] #radial return mapping
                    # print(self.s[i,j].shape, self.s.shape, self.q.shape)
                    self.q[i,j] =  np.sqrt(3./2. * np.sum(self.s[i,j] * self.s[i,j]))
                    self.ε_e[i,j] = self.s[i,j]/(2.*self.G) + ε_e_v_trial[i,j]*self.identity
                    self.Φ[i,j] = 0. #for numerical stability
        print(f'plastičnih točk: {plast} / {self.Ne*self.Nh}')

    def global_loop(self, load_incr=1, tol=None):
        if tol is not None: self.tol = tol
        self.load_incr = load_incr
        self.residuals = np.empty(self.N_iter)
        self.err_elastic = np.empty(self.N_iter)
        self.external_f()
        self.f_int[:] = 0.
        self.r[:] = self.f_int - self.f_ext #force residual
        fnorm = np.linalg.norm(self.f_ext)
        print(fnorm)

        res=self.tol+1
        k=0
        while res > self.tol and k < self.N_iter:
            Δu = spsolve(self.K_glob, self.r)
            self.u += Δu
            # self.u += Δu.reshape(self.vozlisca.shape + (Nps,))
            Δu, theta = Δu.reshape(self.vozlisca.shape).T
            Δε = self.strain_EB(Δu)

            # self.local_loop(Δε)
            self.local1d(Δε) # euler bernoulli

            self.internal_f()
            self.r[:] = self.f_int - self.f_ext
            res = np.linalg.norm(self.r)# / fnorm
            self.residuals[k] = res
            err = np.linalg.norm(self.u.reshape(self.vozlisca.shape).T[0]-self.upogibnica()) / np.linalg.norm(self.upogibnica())
            self.err_elastic[k] = err
            k+=1
            print(f'iteracija {k}, residual: {res:.3g}, {res/fnorm:.3g}, {err:.3g}')

    def incr_load(self, N_incr, max_incr=1):
        for i in range(N_incr):
            load_incr = (i+1)/N_incr * max_incr
            self.global_loop(load_incr=load_incr)
            print(f'end global loop, load increment: {load_incr:.3g}, {i+1} / {N_incr}')
            input('press enter to continue')
    
    def external_f(self):
        self.f_ext[:] = self.load_incr * self.F_glob #[F_glob.reshape((Ne+1,2)) for i in range(self.Nh)]

    def internal_f(self):
        '''izračun notranjih sil iz napetosti za euler bernoulli nosilec'''
        # self.f_int[:] = -self.K_glob @ self.u #elastično - nepopravljeni pomiki!

        self.f_int[:] = 0.
        for i in range(self.Ne): #upogibni moment iz napetosti
            self.M[i] = 2*np.trapz(self.σ[i,:,0,0] * self.y_list, self.y_list) #symmetry

        k, n = linfit(self.integ_tocke[:,0], self.M)
        self.T[:] = k
        # D = findiff.FinDiff(0, hx, 1, acc=4)
        # self.T[:] = D(self.M)
        # self.T[:] = np.gradient(self.M, self.hx, edge_order=1) #prečna sila iz momenta
        # for i in range(1, self.Ne-1): # prečna sila
        #     self.T[i] = (-self.M[i-1] + self.M[i+1]) / (2*self.hx)
        # self.T[0] = (-self.M[0] + self.M[1]) / self.hx
        # self.T[-1] = (-self.M[-2] + self.M[-1]) / self.hx
        # self.T[0] = (-3*self.M[0] + 4*self.M[1] - self.M[2]) / (2*self.hx)
        # self.T[-1] = (self.M[-3] - 4*self.M[-2] + 3*self.M[-1]) / (2*self.hx)
        
        for i in range(self.Ne):
            self.f_int[2*i] += self.T[i]
            self.f_int[2*(i+1)] += self.T[i]
            self.f_int[2*i+1] += self.M[i]
            self.f_int[2*(i+1)+1] += self.M[i]

        # Dirichlt BC
        if self.RP[0] == 1 or self.RP[0] == 2: #prvo vozlisce je fiksno
            self.f_int[0] = self.u[0]
        if self.RP[0] == 2: #prvo vozlisce konzola
            self.f_int[1] = self.u[1]
        if self.RP[1] == 1 or self.RP[1] == 2: #zadnje vozlisce je fiksno
            self.f_int[-2] = self.u[-2]
        if self.RP[1] == 2: #zadnje vozlisce konzola
            self.f_int[-1] = self.u[-1]

    def strain_EB(self, u):
        '''izračun deformacij iz pomikov za euler bernoulli nosilec z metodo končnih razlik'''
        ε = np.zeros_like(self.ε, dtype=float)
        for i in range(1,self.Ne-1):
            kappa = (u[i-1] - 2.*u[i] + u[i+1]) / self.hx**2 #ukrivljenost
            ε[i,:,0,0] = -self.y_list * kappa
        ε[0,:,0,0] = -self.y_list * (2.*u[0] - 5.*u[1] + 4.*u[2] - u[3]) / self.hx**2
        ε[-1,:,0,0] = -self.y_list * (-u[-4] + 4.*u[-3] - 5.*u[-2] + 2.*u[-1]) / self.hx**2
        return ε

    def plot_tez_voz(self, val, label=None, n_mids=10, lw=5, markersize=100):
        '''plottanje težiščnice in vozliščnih vrednosti'''
        fig,ax = plt.subplots(figsize=(16,10))
        cmap = plt.get_cmap('turbo')
        norm = plt.Normalize(vmin=np.min(val), vmax=np.max(val))

        ln=colored_line(self.vozlisca[:,0], self.vozlisca[:,1]+self.U, c=val, ax=ax, cmap=cmap, norm=norm, n_mids=n_mids, linewidths=lw)
        sc=ax.scatter(self.vozlisca[:,0], self.vozlisca[:,1]+self.U, c=val, cmap=cmap, norm=norm, s=markersize)

        fig.colorbar(ln, label=label)
        ax.set_xlabel('$x$ [mm]')
        ax.set_ylabel('$y$ [mm]')
        plt.show()

    def plot_tez_integ(self, val, label=None, n_mids=10, lw=5, markersize=100):
        '''plottanje težiščnice in vrednosti v integracijskih točkah'''
        val_voz =  np.zeros(self.Nv) #integracijske vrednosti lin interpolirane na vozliščih
        for i in range(1,self.Ne): val_voz[i] = (val[i-1] + val[i])/2
        val_voz[0] = val[0] - (val[1]-val[0])/2
        val_voz[-1] = val[-1] + (val[-1]-val[-2])/2

        fig,ax = plt.subplots(figsize=(16,10))
        cmap = plt.get_cmap('turbo')
        norm = plt.Normalize(vmin=np.min(val_voz), vmax=np.max(val_voz))

        ln=colored_line(self.vozlisca[:,0], self.vozlisca[:,1]+self.U, c=val_voz, ax=ax, cmap=cmap, norm=norm, n_mids=n_mids, linewidths=lw)
        sc=ax.scatter(self.vozlisca[:,0], self.vozlisca[:,1]+self.U, c=val_voz, cmap=cmap, norm=norm, s=markersize)

        sc_integ=ax.scatter(self.integ_tocke[:,0], self.integ_tocke[:,1]+self.u_integ, c=val, norm=norm, marker='*', cmap=cmap, s=markersize)
        fig.colorbar(ln, label=label)
        ax.set_xlabel('$x$ [mm]')
        ax.set_ylabel('$y$ [mm]')
        plt.show()

    def plot_cross_voz(self, val, label=None, levels=50):
        '''plottanje celega prereza z vozliščnimi vrednostmi'''
        fig,ax = plt.subplots(figsize=(16,10))
        cmap = plt.get_cmap('turbo')
        norm = plt.Normalize(vmin=np.min(val), vmax=np.max(val))

        cf=ax.contourf(self.X_defor, self.Y_defor, val, levels=levels, cmap=cmap, norm=norm)
        plt.colorbar(cf, label=label)
        ax.set_xlabel('$x$ [mm]')
        ax.set_ylabel('$y$ [mm]')
        plt.show()

    def plot_cross_integ(self, val, label=None, levels=50):
        '''plottanje celega prereza z vrednostmi v integracijskih točkah, interpolirano na vozlišča'''
        val_voz =  np.empty((self.Nv,self.Nh)) #integracijske vrednosti lin interpolirane na vozliščih
        for i in range(1,self.Ne): val_voz[i,:] = (val[i-1,:] + val[i,:])/2
        val_voz[0,:] = val[0,:] - (val[1,:]-val[0,:])/2
        val_voz[-1,:] = val[-1,:] + (val[-1,:]-val[-2,:])/2

        fig,ax = plt.subplots(figsize=(16,10))
        cmap = plt.get_cmap('turbo')
        norm = plt.Normalize(vmin=np.min(val_voz), vmax=np.max(val_voz))

        cf=ax.contourf(self.X_defor, self.Y_defor, val_voz, levels=levels, cmap=cmap, norm=norm)
        plt.colorbar(cf, label=label)
        ax.set_xlabel('$x$ [mm]')
        ax.set_ylabel('$y$ [mm]')
        plt.show()

def remove_neg(a,b):
    # if a[0] <= 0:
    #     a = a[1:]
    #     b = b[1:]
    #     return remove_neg(a,b)
    # else: return a,b
    for i in range(len(a)): #remove ALL negative values
        if a[i] <= 0:
            a = np.delete(a, i)
            b = np.delete(b, i)
            return remove_neg(a,b)
    return a,b
def linfit(x, y):
    '''linearna regresija'''
    A = np.vstack([x, np.ones(len(x))]).T
    k, n = np.linalg.lstsq(A, y, rcond=None)[0]
    return k,n

if __name__ == "__main__":
    # PODATKI
    Nv = 31 # št. vozlišč
    Ne = Nv-1 #št el
    Nh = 10 # št. točk po višini prereza

    # parametri za Newton-Raphson-ovo metodo
    tol = 1e-3
    Nmax = 100

    # # plastičnost - krivulja
    # εn = np.linspace(0,2,1000)
    # σn=[] #inž. napetost, nominalna
    # for i in range(len(εn)):
    #     σ_e = E*εn[i] #elastično
    #     if σ_e <= R_p: σn.append(σ_e)
    #     # else: σn.append(R_p) #idealno plastično -> Et=0
    #     else: σn.append(R_p + Et*(εn[i]-R_p/E)) #konstanten tangentni modul
    # σn = np.asarray(σn)

    # # # podatki iz LNMS
    # # # ROČNO pravi modul elastičnosti E in meja plastičnosti R_p!
    # # εn = np.load('./Modelska_II/zaklj_nal/epsilon.npy')
    # # σn = np.load('./Modelska_II/zaklj_nal/sigma.npy')

    # σt = σn*(1+εn) #Cauchy, prava napetost/true stress
    # ε_log = np.log(1+εn) #logaritemska deformacija
    # ε_pl = ε_log - σt/E #plastična deformacija
    # # ε_pl, σt = remove_neg(ε_pl, σt)

    # plt.plot(εn, σn, 'k-', label='$\sigma_{ing}$')
    # plt.plot(ε_pl, σt, 'b-', label='$\sigma_{true}$')
    # plt.legend()
    # plt.xlabel('$\\varepsilon$ [/]')
    # plt.ylabel('$\\sigma$ [MPa]')
    # plt.show()

    sol = solver(Ne, Nh, Nmax, tol)
    sol.discretization()

    # # reševanje sistema enačb
    # t2=timer()
    # u_elast = spsolve(K_glob, F_glob)
    # t3=timer()
    # U_elast, Θ_elast = u_elast.reshape(sol.vozlisca.shape).T
    # print(f'časi reševanja: sparse {(t3-t2):.3g} s')
    # print(f'maksimalni poves: {np.max(np.abs(U_elast)):.5g} mm')

    # plastično
    #############################################################################

    load = 1e-1
    sol.N_iter = 100
    sol.tol = 1e-1

    u_elast = spsolve(sol.K_glob, load*sol.F_glob)
    U_elast, Θ_elast = u_elast.reshape(sol.vozlisca.shape).T

    sol.global_loop(load_incr=load)
    # sol.incr_load(N_incr=10, max_incr=load)

    U, sol.Θ = sol.u.reshape(sol.vozlisca.shape).T
    sol.U = -U
    # sol.U,sol.Θ = U_elast, Θ_elast
    sol.u_integ = np.array([(sol.U[i] + sol.U[i+1])/2 for i in range(sol.Ne)]) #poves lin interpoliran na sredini KE

    U_all = np.repeat(sol.U,sol.Nh).reshape(sol.Nv,sol.Nh) #poves razširjen na cel prerez
    sol.X_defor, sol.Y_defor = sol.X_defor, sol.Y_defor + U_all
    
    u_ex = -sol.upogibnica()
    print('elastična teorija:')
    print(f'sigma_max = {load*sol.F0*sol.L*sol.h/2/sol.I:.3g}, izračun: {np.max(np.abs(sol.q)):.3g} MPa')
    print(f'u_max     = {load*sol.F0*sol.L**3/3/sol.E/sol.I:.3g}, izračun: {np.max(sol.u):.3g} mm')
    print(f'M_max     = {load*sol.F0*sol.L:.3g}, izračun: M_max = {np.max(sol.M):.3g} mm')
    print(f'T_max     = {load*sol.F0:.3g}, izračun: {np.max(sol.T):.3g} mm')
    print(f'napaka    : {np.linalg.norm(u_ex-U)/np.linalg.norm(u_ex):.3g}, {np.linalg.norm(U-U_elast):.3g}')
    # input('press enter to continue')

    plt.plot(sol.M, 'x-',label='$M$')
    plt.plot(sol.T, 'x-', label='$T$')
    plt.legend()
    plt.show()
    
    iters = np.arange(sol.N_iter)+1
    plt.loglog(iters,sol.residuals, 'kx-')
    plt.loglog(iters,sol.err_elastic, 'bx-')
    plt.show()

    # plot_tez_integ(sol.M, label='$M$')
    # plot_tez_integ(sol.T, label='$T$')
    sol.plot_cross_voz(np.abs(U_all), label='$|U|$')
    sol.plot_cross_integ(sol.Φ, label='$\Phi$')
    sol.plot_cross_integ(sol.q, label='$\sigma_p$')
    # sol.plot_cross_integ(sol.σ[:,:,0,0], label='$\sigma_{xx}$')
    sol.plot_cross_integ(sol.ε_p_eq, label='$\\varepsilon_p^{eq}$')

    sol.plot_cross_integ(sol.ε[:,:,0,0], label='$\\varepsilon_{xx}$')

    tocka = 0 #0 nevtralana os ali -1 za prikaz ektremnih vrednosti na nosilcu - upoštevamo simetrijo
    # tocka = 0 #N//2 za nevtralana os ali 0 oz. -1 za prikaz ektremnih vrednosti na nosilcu - cel prerez

    sol.plot_tez_voz(np.abs(U), label='$|u|$')
    # sol.plot_tez_integ(np.abs(u_integ), label='$|u_{integ}|$')
    sol.plot_tez_integ(sol.q[:,tocka], label='$\sigma_m(N.O.)$')

    tocka=-1
    sol.plot_tez_integ(sol.q[:,tocka], label='$\sigma_m(y_{max})$')
    sol.plot_tez_integ(sol.Φ[:,tocka], label='$\Phi$')