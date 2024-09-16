# coding: utf-8
#SIMPLE_demo, a demonstration of fluid flow solver using SIMPLE algorithm 
#Copyright (C) 2020  Boštjan Mavrič

from __future__ import division
import matplotlib.pyplot as plt
import math as mt
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sp_linalg
import scipy.integrate as integrate
from IPython.display import display, clear_output
from matplotlib.patches import FancyArrowPatch
from matplotlib import rc
from timeit import default_timer as timer

rc({"text.usetex": True,
    'family': 'serif', 'serif': ['Computer Modern']})
def pause():
    input("Press the <ENTER> key to continue...")

# Preslikava med indeksi premaknjene mreže kontrolnih volumnov na dolge vektorje za reševanje sistema linearnih enačb.
# / Mapping of the indices from staggered grid to the "long" vectors used to solve systems of linear equations.
def glob_ind_u(i, j, N):
    return (i) + j*(N+1)

def glob_ind_v(i, j, N):
    return (i) + j*N

def glob_ind_p(i, j, N):
    return (i) + j*N

# Interpolacija hitrostnih vrednosti na položaje tlačnih točk.
# / Interpolate velocity values to the positions of pressure variables.
def interpolate_to_pressure_pos(us, vs, vec, N):
    for i in range(0, N):
        for j in range(0, N):
            vec[0, i, j]=(us[i, j]+ us[i+1, j])/2
            vec[1, i, j]=(vs[i, j]+ vs[i, j+1])/2

# Preslikaj  vrednosti iz dolgih vektorjev na matrike, ki ustrezajo zamaknjeni mreži kontrolnih volumov. /
# / Map values from "long" vectors to the arrays that correspond to staggered grid.
def reshuffle_us(u_sol, us, N):
    for i in range(0, N+1):
        for j in range(0, N):
            us[i, j]=u_sol[glob_ind_u(i, j, N)]

def reshuffle_vs(v_sol, vs, N):
    for i in range(0, N):
        for j in range(0, N+1):
            vs[i, j]=v_sol[glob_ind_v(i, j, N)]

def reshuffle_p(p_sol, ps, N):
    for i in range(0, N):
        for j in range(0, N):
            ps[i, j]=p_sol[glob_ind_p(i, j, N)]

# Funkcija, ki reši sistem enačb.
# / Main function used to solve the system of equations.
def SIMPLE(N, FD_kind, Re, alpha_p, alpha_vel, Nmax_iters=200,v_lid=1,rho=1,tol=1e-3,b=1):
    t_start = timer()

    # monitoring residualov
    p_corr_norm=np.zeros(Nmax_iters)
    u_corr_norm=np.zeros(Nmax_iters)
    v_corr_norm=np.zeros(Nmax_iters)
    continuity=np.zeros(Nmax_iters)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18,6))
    ax1.set_xlabel('$x$')
    ax1.set_ylabel('$y$')
    ax2.set_xlabel('$N_{iter}$')
    ax2.set_ylabel('scaled residual')
    ax3.set_xlabel('$x$, $y$')
    ax3.set_ylabel('$u$, $v$')
    plt.ion()
    plt.show()
    
    # Shranjevanje spremenljivk.
    # / Storage for unknowns.
    ############################

    # Shramba za u. 
    # / Storage of u.
    # u ima eno vrstico spremenljivk več vzdolž x kot vzdolž y.
    # / For u, there is one more row of variables along x axis than along y.
    Nall_u=(N)*(N+1) #število neznank
    us=np.zeros([N+1, N]) #matrika rešitev
    us_old=np.zeros([N+1, N]) #prejšnja rešitev
    ds_u=np.zeros([N+1, N]) #diskretizirana enačba - inverz diagonale za tlačni popravek (L⁻¹)
    rhs_u=np.zeros(Nall_u) #desna stran N-S en. - diskretizacija tlaka (vektor G)
    u_sol=np.zeros(Nall_u) #vektor neznank, vrne solver -> transformiramo v matriko rešitev

    # Shramaba za v
    # / Storage of v.
    # v ima eno eno vrstico spremenljivk vzdolž y osi kot vzdolž x.
    # / There is one more row of v variables along y axis than along x.
    Nall_v=(N+1)*(N) #=Nall_u
    vs=np.zeros([N, N+1])
    vs_old=np.zeros([N, N+1])
    ds_v=np.zeros([N, N+1])
    rhs_v=np.zeros(Nall_v)
    v_sol=np.zeros(Nall_v)

    # Shramba za p.
    # / Storage of p.
    Nall_p=(N)*(N) #zamaknjena mreža -> za N manj neznank
    ps=np.zeros([N, N]) # trenutna rešitev
    ps_prime=np.zeros([N, N]) #tlačni popravek v matrični obliki
    vel_at_p=np.zeros([2, N, N]) #interpolirane hitrosti na mesta izračuna tlaka (zamaknjena mreža)
    vel_norm_at_p=np.zeros([N, N]) #normirana hitrost -> velikost vektorja hitrosti?

    # Za stabilizacijo tlačnega popravka potrebujemo eno enačbo več kot je tlačnih vrednosti.
    # / To stabilize pressure correction equation we need on more unknown for regularization of the system.
    rhs_p=np.zeros(Nall_p+1) #diskretizacija divergence hitrosti - desna stran tlačnega popravka   
    p_sol=np.zeros(Nall_p+1) #vektor neznank - tlačni popravki
    p_old=np.zeros(Nall_p+1) #prejšnji popravek

    # Shranjevanje matrik.
    # / Storage for matrices.
    # Redke matrike sestavimo v obliki "seznama seznamov". Preden jih uporabimo v reševalniku jih spremenimo v format CSR.
    # / We use sparse matrices in "List-of-lists" for construction. Before using them in the solver, they are converted to CSR format.
    # globalna matrika H sestoji iz H1 in H2 -> SIMPLE uporablja le diagonalni del L=diag(H), ki je shranjen pod rhs_u [diag(H1)] in rhs_v [diag(H2)]
    sparse_u=sp.lil_matrix((Nall_u, Nall_u)) #H1
    sparse_v=sp.lil_matrix((Nall_v, Nall_v)) #H2
    sparse_p=sp.lil_matrix((Nall_p+1, Nall_p+1)) #globalna matrika sistema enačb tlačnega popravka - D (L⁻¹ G p')

    # Generirajmo zamaknjeno mrežo. - KOORDINATE RAČUNSKIH TOČK
    # / Generation of staggered grid.
    ###############################
    h=1./N
    A=h*b #presek KV za 2D problem, b-debelina KV

    # Točke z vrednostjo u.
    # / Points with known value of u.
    points_u=np.zeros([2, Nall_u])
    for i in range(0, N+1):
        for j in range(0, N):
            Ip=glob_ind_u(i, j, N)
            points_u[0, Ip]=h*i
            points_u[1, Ip]=h/2 + h*j

    # Točke z vrednostjo v.
    # / Points with known value of v.
    points_v=np.zeros([2, Nall_v])
    for i in range(0, N):
        for j in range(0, N+1):
            Ip=glob_ind_v(i, j, N)
            points_v[0, Ip]=h/2 + h*i
            points_v[1, Ip]=h*j

    # Točke z vrednostjo tlaka.        
    # / Points with known value of p.        
    points_p=np.zeros([2, Nall_p])
    points_p_contour=np.zeros(N)
    for i in range(0, N):
        points_p_contour[i]=h/2 + h*i
        for j in range(0, N):
            Ip=glob_ind_p(i, j, N)
            points_p[0, Ip]=h/2 + h*i
            points_p[1, Ip]=h/2 + h*j
    
    
    # Izraze poenostavimo za konstantne vrednosti gostote in transportnih koeficientov.
    # / The expressions are simplified by assuming constant density and constant transport coefficients.
    Gamma=1e0/Re #brezdimenzijska oblika N-S en
   
    # Difuzijski koeficienti iz V&M učbenika.
    # / Diffusion coefficents from V&M book.
    Dw=Gamma/h
    De=Gamma/h
    Ds=Gamma/h
    Dn=Gamma/h
    
    # Začetek iteracijske zanke.
    # / Start of the the iteration loop.
    ############################
    itr=0
    residual = 1
    while itr < Nmax_iters and residual > tol:
        Pe_u=0.
        
        # Pripravi sistem enačb za poskusno hitrost (u, v).
        # / Set up the system of equations for guess velocity (u, v).        
        ###################################
        # u-komponenta hitrosti, notranjosti.
        # / u-component of velocity, interior.
        for i in range(1, N):
            for j in range(1, N-1):
                # Določi potrebne indekse točk.
                # / Determine indices of points under consideration.
                Ip=glob_ind_u(i, j, N)
                Ie=glob_ind_u(i+1, j, N)
                Iw=glob_ind_u(i-1, j, N)
                In=glob_ind_u(i, j+1, N)
                Is=glob_ind_u(i, j-1, N)

                # Določi advektivne tokove.
                # / Determine advective fluxes. -> central diferences
                Fw=rho*(us[i, j] + us[i-1, j])/2
                Fe=rho*(us[i+1, j] + us[i, j])/2
            
                Fs=rho*(vs[i, j] + vs[i-1, j])/2
                Fn=rho*(vs[i, j+1] + vs[i-1, j+1])/2

                # Določi koeficiente v linearnem sistemu enačb.
                # / Determine the coefficents used in linear systems of equations.
                if(FD_kind=="central"):
                    a_w=Dw+Fw/2
                    a_e=De-Fe/2
                    a_s=Ds+Fs/2
                    a_n=Dn-Fn/2
                elif(FD_kind=="upwind"):
                    a_w=Dw+max([Fw, 0.])
                    a_e=De+max([-Fe, 0.])
                    a_s=Ds+max([Fs, 0.])
                    a_n=Dn+max([-Fn, 0.])
                elif(FD_kind=="hybrid"):
                    a_w=max([Fw, (Dw+Fw/2), 0])
                    a_e=max([-Fe, (De-Fe/2), 0])
                    a_s=max([Fs, (Ds+Fs/2), 0])
                    a_n=max([-Fn, (Dn-Fn/2), 0])

                a_p=a_w+a_e+a_s+a_n +(Fe-Fw)+(Fn-Fs)

                # Nadzorujmo Pecletovo število.
                # / Monitor the Peclet number.
                Pe_u=max([Pe_u, max([Fw, Fe, Fs, Fn])/Dw])
                
                # Nastavimo matrične elemente.
                # / Populate the matrix elements.
                sparse_u[Ip, In]=-a_n
                sparse_u[Ip, Iw]=-a_w
                sparse_u[Ip, Ip]=a_p
                sparse_u[Ip, Ie]=-a_e
                sparse_u[Ip, Is]=-a_s

                # Shranimo vrednosti za tlačni popravek.
                # / Store values for pressure correction.
                ds_u[i, j]=h/a_p
            
                # Nastavimo vrednost desne strani za u komponento.
                # / Set up RHS for u-component.
                rhs_u[Ip]=A*(ps[i-1, j] - ps[i, j])

        # u-kompnenta hitrosti v točkah v stiku s severno/južno steno.
        # / u-component of velocity in cells in contact with north/south walls.
        for i in range(1, N):
            for j in [0, N-1]:
                
                Ip=glob_ind_u(i, j, N)
                Ie=glob_ind_u(i+1, j, N)
                Iw=glob_ind_u(i-1, j, N)
                In=glob_ind_u(i, j+1, N)
                Is=glob_ind_u(i, j-1, N)
                
                Fw=rho*(us[i, j] + us[i-1, j])/2
                Fe=rho*(us[i+1, j] + us[i, j])/2
            
                Fs=rho*(vs[i, j] + vs[i-1, j])/2
                Fn=rho*(vs[i, j+1] + vs[i-1, j+1])/2

                if(FD_kind=="central"):
                    a_w=Dw+Fw/2
                    a_e=De-Fe/2
                    a_s=Ds+Fs/2
                    a_n=Dn-Fn/2
                elif(FD_kind=="upwind"):
                    a_w=Dw+max([Fw, 0.])
                    a_e=De+max([-Fe, 0.])
                    a_s=Ds+max([Fs, 0.])
                    a_n=Dn+max([-Fn, 0.])
                elif(FD_kind=="hybrid"):
                    a_w=max([Fw, (Dw+Fw/2), 0])
                    a_e=max([-Fe, (De-Fe/2), 0])
                    a_s=max([Fs, (Ds+Fs/2), 0])
                    a_n=max([-Fn, (Dn-Fn/2), 0])

                sparse_u[Ip, Iw]=-a_w
                sparse_u[Ip, Ie]=-a_e

                # Stena pri y=0.
                # / Wall at y=0.
                if(j==0):
                    a_s=a_s+Ds
                    a_p=a_w+a_e+a_s+a_n +(Fe-Fw)+(Fn-Fs)
                    sparse_u[Ip, Ip]=a_p
                    sparse_u[Ip, In]=-a_n #a_s = 0

                    # Nastavimo vrednost v vektorju desnih strani.
                    # / Set up RHS for u-component.
                    rhs_u[Ip]=A*(ps[i-1, j] - ps[i, j])
                    
                # Stena pri y=1, kjer je gibajoč pokrov.
                # / Wall at y=1, the position of moving lid.
                elif(j==N-1):
                    a_n=a_n+Dn
                    a_p=a_w+a_e+a_s+a_n + (Fe-Fw)+(Fn-Fs)
                    sparse_u[Ip, Ip]=a_p
                    sparse_u[Ip, Is]=-a_s #a_n = 0

                    # Nastavimo desne strani za u-komponento. Dodatni člen povzroči gibanje tekočine.
                    # / Set up RHS for u-component, the last term provides the driving.
                    rhs_u[Ip]=A*(ps[i-1, j] - ps[i, j]) + a_n*v_lid 
                
                ds_u[i, j]=h/a_p               
            
        # u kompnenta hitrosti v celicah na stiku z vzhodno/zahodno steno.
        # / u-component of velocity in cells in contact with east/west walls.
        for i in [0, N]:
            for j in range(0, N):

                Ip=glob_ind_u(i, j, N)

                # Dirichletovi robni pogoji.
                # / Dirichlet boundary conditions.
                sparse_u[Ip, Ip]=1
            
                ds_u[i, j]=0
            
                # Nastavimo vrednosti v vektorju desnih strani za u-komponento.
                # / Set up RHS for u-component.
                rhs_u[Ip]=0

        # Ponovimo postopek še za v-komponento.
        # / The same procedure for v-component.
        ###################################
        # v-komponenta hitrosti , notranjost.
        # v-component of velocity, interior.
        for i in range(1, N-1):
            for j in range(1, N):
                
                Ip=glob_ind_v(i, j, N)
                Ie=glob_ind_v(i+1, j, N)
                Iw=glob_ind_v(i-1, j, N)
                In=glob_ind_v(i, j+1, N)
                Is=glob_ind_v(i, j-1, N)
            
                Fw=rho*(us[i, j] + us[i, j-1])/2
                Fe=rho*(us[i+1, j] + us[i+1, j-1])/2
            
                Fs=rho*(vs[i, j-1] + vs[i, j])/2
                Fn=rho*(vs[i, j] + vs[i, j+1])/2

                if(FD_kind=="central"):
                    a_w=Dw+Fw/2
                    a_e=De-Fe/2
                    a_s=Ds+Fs/2
                    a_n=Dn-Fn/2
                elif(FD_kind=="upwind"):
                    a_w=Dw+max([Fw, 0.])
                    a_e=De+max([-Fe, 0.])
                    a_s=Ds+max([Fs, 0.])
                    a_n=Dn+max([-Fn, 0.])
                elif(FD_kind=="hybrid"):
                    a_w=max([Fw, (Dw+Fw/2), 0.])
                    a_e=max([-Fe, (De-Fe/2), 0.])
                    a_s=max([Fs, (Ds+Fs/2), 0.])
                    a_n=max([-Fn, (Dn-Fn/2), 0.])
                
                a_p=a_w+a_e+a_s+a_n +(Fe-Fw) +(Fn-Fs)

                sparse_v[Ip, In]=-a_n
                sparse_v[Ip, Iw]=-a_w
                sparse_v[Ip, Ip]=a_p
                sparse_v[Ip, Ie]=-a_e
                sparse_v[Ip, Is]=-a_s

                ds_v[i, j]=h/a_p

                # Nastavimo vrednost v vektorju desnih strani za v-komponento.
                # / Set up RHS for v-component.
                rhs_v[Ip]=A*(ps[i, j-1] - ps[i, j])

        # v-komponenta hitrosti v celicah v stiku z vzhodno/zahodno steno.
        # / v-component of velocity in cells in contact with east/west walls.
        for i in [0, N-1]:
            for j in range(1, N):
            
                Ip=glob_ind_v(i, j, N)
                Ie=glob_ind_v(i+1, j, N)
                Iw=glob_ind_v(i-1, j, N)
                In=glob_ind_v(i, j+1, N)
                Is=glob_ind_v(i, j-1, N)
            
                Fw=rho*(us[i, j] + us[i-1, j])/2
                Fe=rho*(us[i+1, j] + us[i, j])/2
            
                Fs=rho*(vs[i, j] + vs[i-1, j])/2
                Fn=rho*(vs[i, j+1] + vs[i-1, j+1])/2

                if(FD_kind=="central"):
                    a_w=Dw+Fw/2
                    a_e=De-Fe/2
                    a_s=Ds+Fs/2
                    a_n=Dn-Fn/2
                elif(FD_kind=="upwind"):
                    a_w=Dw+max([Fw, 0.])
                    a_e=De+max([-Fe, 0.])
                    a_s=Ds+max([Fs, 0.])
                    a_n=Dn+max([-Fn, 0.])
                elif(FD_kind=="hybrid"):
                    a_w=max([Fw, (Dw+Fw/2), 0])
                    a_e=max([-Fe, (De-Fe/2), 0])
                    a_s=max([Fs, (Ds+Fs/2), 0])
                    a_n=max([-Fn, (Dn-Fn/2), 0])

                sparse_v[Ip, In]=-a_n
                sparse_v[Ip, Is]=-a_s
            
                if(i==0):
                    a_w=a_w + Dw
                    a_p=a_w+a_e+a_s+a_n +(Fe-Fw)+(Fn-Fs)
                    a_p=a_p 
                    sparse_v[Ip, Ip]=a_p
                    sparse_v[Ip, Ie]=-a_e
                elif(i==N-1):
                    a_e=a_e + De
                    a_p=a_w+a_e+a_s+a_n +(Fe-Fw)+(Fn-Fs)
                    a_p=a_p
                    sparse_v[Ip, Ip]=a_p
                    sparse_v[Ip, Iw]=-a_w
                
                ds_v[i, j]=h/a_p
            
                # Nastavi vrednosti v vektorju desnih strani za v-komponento.
                # / Set up RHS for v-component.
                rhs_v[Ip]=A*(ps[i, j-1] - ps[i, j])

        # v-komponenta hitrosti v celicah v stiku s severno/južno steno.
        # / v-component of velocity in cells in contact with north/south walls.
        for i in range(0, N):
            for j in [0, N]:

                Ip=glob_ind_v(i, j, N)
            
                sparse_v[Ip, Ip]=1
            
                ds_v[i, j]=0

                # Nastavi vrednosti v vektorju desnih strani za u-komponento.
                # / Set up RHS for u-component.
                rhs_v[Ip]=0
        
        # Rešimo sistema enačb za poskusno hitrost.
        ###################################
        # Rešimo sistem enačb za u.
        # / Solve the system of equations for u.
        u_sol=sp_linalg.spsolve(sparse_u.tocsr(), rhs_u)
        reshuffle_us(u_sol, us, N)

        # Rešimo sistem enačb za v.
        # / Solve the system of equations for v.
        v_sol=sp_linalg.spsolve(sparse_v.tocsr(), rhs_v)
        reshuffle_vs(v_sol, vs, N)
        
        # Izračunamo tlačni popravek iz poskusne hitrosti.
        ###################################
        # Sestavimo tlačno matriko.
        # / Assemble the pressure matrix.
        for i in range(0, N):
            for j in range(0, N):
                Ip=glob_ind_p(i, j, N)
                Ie=glob_ind_p(i+1, j, N)
                Iw=glob_ind_p(i-1, j, N)
                In=glob_ind_p(i, j+1, N)
                Is=glob_ind_p(i, j-1, N)

                a_e=ds_u[i+1, j]/h
                a_w=ds_u[i, j]/h
                a_n=ds_v[i, j+1]/h
                a_s=ds_v[i, j]/h

                # Obravnavamo robne pogoje.
                # / Treat the boundary conditions.
                if(i==0) :
                    a_w=2*a_w
                else:
                    sparse_p[Ip, Iw]=-a_w

                if(i==N-1):    
                    a_e=2*a_e
                else:
                    sparse_p[Ip, Ie]=-a_e
                
                if(j==0):
                    a_s=2*a_s
                else:
                    sparse_p[Ip, Is]=-a_s
                
                if(j==N-1):
                    a_n=2*a_n
                else:
                    sparse_p[Ip, In]=-a_n
                
                a_p=a_e+a_w+a_n+a_s
                sparse_p[Ip, Ip]=a_p
            
                # Izračunajmo vrednosti v vektorju desnih strani.
                # / Determine values in RHS vector for pressure.
                rhs_p[Ip]=h*(us[i,j] - us[i+1,j] + vs[i,j] - vs[i,j+1])   

        # Dodamo regularizacijsko enačbo.
        # / Add the regularization equation sum(p)=0.
        rhs_p[Nall_p]=0
        sparse_p[0:Nall_p, Nall_p]=1e0
        sparse_p[Nall_p, 0:Nall_p]=1e0

        # Rešimo sistem enačb.
        # / Solve the system of equations.
        p_sol=sp_linalg.spsolve(sparse_p.tocsr(), rhs_p)
        reshuffle_p(p_sol, ps_prime, N) #transformiramo dolg vektor v matriko
        
        # Izračunamo hitrostni popravek. -> direkten izračun nove hitrosti!
        ###################################
        # Popravimo vrednosti u.
        # / Correct u values.
        for i in range(1, N):
            for j in range(0, N):
                us[i, j]=us[i, j] + ds_u[i, j]*(ps_prime[i-1, j]-ps_prime[i, j])
      
        # Popravimo vrednosti v.
        # / Correct v values.
        for i in range(0, N):
            for j in range(1, N):
                vs[i, j]=vs[i, j] + ds_v[i, j]*(ps_prime[i, j-1]-ps_prime[i, j])
        
        # Popravimo tlak in hitrost.
        ###################################
        # Posodobimo vrednosti tlaka.
        # / Update pressure.
        ps=ps+alpha_p*ps_prime

        # Podrelaksacija popravkov hitrosti.
        # / Under-relax velocity corrections.
        us=alpha_vel*us + (1.-alpha_vel)*us_old
        vs=alpha_vel*vs + (1.-alpha_vel)*vs_old
        
        # Izračun residualov
        ########################################
        # popravki, skalirani na prejšnjo vrednost
        p_corr_norm[itr]=np.linalg.norm(ps_prime)/(np.linalg.norm(ps))
        u_corr_norm[itr]=np.linalg.norm(us-us_old)/(np.linalg.norm(us_old))
        v_corr_norm[itr]=np.linalg.norm(vs-vs_old)/(np.linalg.norm(vs_old))

        # sedaj lahko shranimo shranimo stare vrednosti
        p_old=ps.copy() #to je popravek!
        us_old=us.copy()
        vs_old=vs.copy()
        
        # kontinuiteta - absolutni pogrešek (res == 0)
        residual=0
        for i in range(0,N):
            for j in range(0, N):
                residual += abs(-us[i,j] + us[i+1,j] - vs[i,j] + vs[i,j+1])*h
        continuity[itr] = residual

        # print(f'iteracija {itr+1}/{Nmax_iters}:')
        # print(f'\t p = {ps[N//2,N//2]:.5f}')
        # print(f'\t u = {us[N//2,N//2]:.5f}')
        # print(f'\t v = {vs[N//2,N//2]:.5f}')

        # Priprava podatkov za izris.
        # / Extract data for plotting.
        interpolate_to_pressure_pos(us, vs, vel_at_p, N)
        for i in range(0, N):
            for j in range(0, N):
                vel_norm_at_p[i, j]=np.linalg.norm(vel_at_p[:, i, j])

        pmin=np.amin(ps)
        pmax=np.amax(ps)

        ax1.cla()
        ax2.cla()
        ax3.cla()

        ax1.set_xlabel('$x$')
        ax1.set_ylabel('$y$')
        ax2.set_xlabel('$N_{iter}$')
        ax2.set_ylabel('scaled residual')
        ax3.set_xlabel('$x$, $y$')
        ax3.set_ylabel('$u$, $v$')

        ax1.set_title("Streamlines and pressure field")
        
        ax1.contourf(points_p_contour, points_p_contour, ps.transpose(), levels=np.linspace(pmin, pmax, 30), cmap='YlGnBu_r')
        ax1.streamplot(points_p[0, 0:N], points_p[1, 0:Nall_p:N], vel_at_p[0, :, :].transpose(), vel_at_p[1, :, :].transpose(), color=vel_norm_at_p.transpose(), density=0.5, linewidth=3., cmap='hot', broken_streamlines=True) #hitreje če so prekinjene tokovnice
        ax1.axhline(0.5,c='C0')
        ax1.axvline(0.5,c='C1')

        ax1.set_aspect('equal')
        ax1.set_xlim(0, 1)
        ax1.set_ylim(0, 1)
        
        ax2.set_title("Relative magnitude of corrections")
        ax2.set_yscale('log')
        ax2.plot(p_corr_norm[1:itr], 'g',label="Pressure")
        ax2.plot(u_corr_norm[1:itr], 'b',label="Velocity x")
        ax2.plot(v_corr_norm[1:itr], 'c',label="Velocity y")
        # ax2.plot(continuity[1:itr]/continuity[itr-1], label='Continuity')
        ax2.plot(continuity[1:itr]/continuity[0], 'k',label='Continuity')
        ax2.legend()
        ax2.grid()

        ax3.plot(points_p_contour, vs[:, int(N/2)+1], label="v @ y=0.5")
        ax3.plot(points_p_contour, us[int(N/2)+1, :], label="u @ x=0.5")
        ax3.set_title("Velocity profiles")
        ax3.legend()
        ax3.grid()
        #plt.draw()
        fig.suptitle(f'iteration {itr+1}/{Nmax_iters}')
        clear_output(wait=True)
        #display(fig)
        fig.show()
        plt.pause(.0001)

        # print(sparse_p.toarray())
        # print(sparse_u.toarray())
        # print(sparse_v.toarray())
        # print(rhs_u)
        # print(rhs_v)

        if residual <= tol: print(f'Solution converged in {itr+1} iterations.')
        elif itr == Nmax_iters-1: print('Number of max iterations reached without convergence.')
        else:
            print('iteration   \t continuity')
            print(f'{itr+1:9d} : \t {residual:.5E}')
        itr+=1

    # Konec iteracij.
    #############################
    plt.ioff()
    plt.show()
    t_end = timer()
    time = t_end-t_start
    print(f'čas simulacije: {time: 10.5f} s')

    # Post proces: izračun sile na pokrov
    F_y_lid = integrate.simpson(ps[:,N-1], dx=h)
    F_y_bot = integrate.simpson(ps[:,0], dx=h)
    F_x_l = integrate.simpson(ps[0,:], dx=h)
    F_x_d = integrate.simpson(ps[N-1,:], dx=h)
    F = np.array([F_y_lid, F_y_bot, F_x_l, F_x_d])
    print(f'sila na pokrov: {F_y_lid: 10.5f} N')
    print(f'sila na dno   : {F_y_bot: 10.5f} N')
    print(f'sila na levo  : {F_x_l: 10.5f} N')
    print(f'sila na desno : {F_x_d: 10.5f} N')

    # Izriši končno rešitev.
    # / Plot final solution.
    fig,ax=plt.subplots(figsize=(8,6.5))
    plt.axis('equal')
    plt.title("Velocity streamlines")
    st=ax.streamplot(points_p[0, 0:N], points_p[1, 0:Nall_p:N], vel_at_p[0, :, :].transpose(), vel_at_p[1, :, :].transpose(), color=vel_norm_at_p.transpose(), density=3, broken_streamlines=False, cmap='turbo')
    for art in ax.get_children(): #remove arrows
        if not isinstance(art, FancyArrowPatch):
            continue
        art.remove()        # Method 1
        #art.set_alpha(0)  # Method 2
    plt.colorbar(st.lines,ax=ax,label='magnitude')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.show()

    fig,ax=plt.subplots(figsize=(8,8))
    plt.axis('equal')
    plt.title("Velocity field") #kaj pa izris v točkah izračuna?
    st=ax.quiver(points_p[0, 0:N], points_p[1, 0:Nall_p:N], vel_at_p[0, :, :].transpose(), vel_at_p[1, :, :].transpose(), scale=10)
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.show()
    
    pmin=np.amin(ps)
    pmax=np.amax(ps)
    
    fig,ax=plt.subplots(figsize=(8,6.5))
    plt.axis('equal')
    plt.title("Pressure field")
    # cf=plt.contourf(points_p_contour, points_p_contour, ps.transpose(), levels=np.linspace(pmin, pmax, 20))
    cf=plt.contourf(points_p_contour, points_p_contour, ps.transpose(), levels=30)
    plt.colorbar(cf,ax=ax)
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.show()

    # fig,ax=plt.subplots(figsize=(8,6.5))
    # plt.axis('equal')
    # plt.title("Pressure field")
    # cf=plt.contourf(points_p[0,:].reshape((N,N)), points_p[1,:].reshape((N,N)), ps.transpose()) #ISTO
    # plt.colorbar(cf,ax=ax)
    # plt.xlabel('$x$')
    # plt.ylabel('$y$')
    # plt.show()

    return F, time, continuity, itr