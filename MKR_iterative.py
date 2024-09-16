# functions for generating and iterative solving system of equations with finite difference method FDM/MKR

import numpy as np

#system generation
def A_plosca(N):
    Ny=Nx=N
    r=1 #polmer cevi - x
    L=1 #dolžina cevi - y
    dx=L/(Nx-1)
    dy=r/(Ny-1)

    # DE, RP
    A = np.zeros(((Nx)*(Ny), (Nx)*(Ny)))
    b = np.zeros(((Nx)*Ny))

    # DE
    for j in range(1,Ny-1):
        for i in range(1,Nx-1):
            ij_index = j*(Nx) + i #- Nx-1
            ij = 2*(1/dx**2 + 1/dy**2)
            ij_dx = -1/dx**2 #- 1/(2*i*dx**2) #polarni dodatek!
            ij_dy = -1/dy**2

            #center
            A[ij_index, ij_index] += ij
            # okolica
            A[ij_index, ij_index-1] += ij_dx #levo
            A[ij_index, ij_index+1] += ij_dx #desno
            A[ij_index, ij_index-Nx] += ij_dy #dol
            A[ij_index, ij_index+Nx] += ij_dy #gor
            
            b[ij_index] = 0 #ni izvirov

    # RP Dirichlet
    T_rp1 = 1
    T_rp2 = 0
    for i in range(Nx):
        #spodaj
        #A[i, i] = 1
        #b[i] = T_rp1
        #zgoraj
        A[-Nx+i, -Nx+i] = 1
        b[-Nx+i] = T_rp1
    for j in range(Ny):
        #levo
        #A[j*Nx, j*Nx] = 1
        #desno
        A[j*Nx+Nx-1, j*Nx+Nx-1] = 1
        if j*dy < r/2: #spodaj
            #b[j*Nx] = T_rp1 #L
            b[j*Nx+Nx-1] = T_rp1 #R
        else: #zgoraj
            #b[j*Nx] = T_rp2 #L
            b[j*Nx+Nx-1] = T_rp2 #R
    # Neumann RP
    for i in range(Nx): #spodaj
        #A[i, i] = 1/(2*dy) #centralna diferenčna shema okoli 1.vrste - vzame 0. in 2. vrsto
        #A[i, i+2*Nx] = -1/(2*dy)
        A[i, i] = 1/(dy) #diferenčna shema naprej okoli 0.vrste - vzame 0. in 1. vrsto
        A[i, i+Nx] = -1/(dy)
        b[i] = 2 #+tok noter, -ven
    for j in range(1,Ny): #levo SYM RP - izolirano
        #A[j*Nx, j*Nx] = 1/(2*dx) #centralna diferenčna shema okoli 1.solpca - vzame 0. in 2. stolpec
        #A[j*Nx, j*Nx+2] = -1/(2*dx)
        A[j*Nx, j*Nx] = 1/(dx) #diferenčna shema naprej okoli 0.stolpca - vzame 0. in 1. stolpec
        A[j*Nx, j*Nx+1] = -1/(dx)
        b[j*Nx] = 0
    # for j in range(0,Ny): #desno (zgoraj) toplotno izolirano, NAMESTO T_rp
    #     #A[j*Nx+Nx-1, j*Nx+Nx-2] = 1/(2*dy) #centralna diferenčna shema okoli zadnjga stolpca - vzame -1. in -3. stolpec
    #     #A[j*Nx+Nx-1, j*Nx+Nx-1] = -1/(2*dy)
    #     A[j*Nx+Nx-1, j*Nx+Nx-2] = -1/(dy) #diferenčna shema nazaj okoli zadnjega stolpca - vzame -1. in -2. stolpec
    #     A[j*Nx+Nx-1, j*Nx+Nx-1] = 1/(dy)
    #     b[j*Nx+Nx-1] = 0
    return A,b

def RP(A,b,insulated,N):
    Nx=Ny=N
    dx=dy=1./(Nx-1)

    # RP Dirichlet
    T_rp1 = 0
    T_rp2 = 1
    for i in range(Nx):
        #spodaj
        #A[i, i] = 1
        #b[i] = T_rp1
        #zgoraj
        A[-Nx+i, -Nx+i] = 1
        b[-Nx+i] = T_rp1
    for j in range(Ny//2):
            #desno spodaj
            A[j*Nx+Nx-1, j*Nx+Nx-1] = 1
            b[j*Nx+Nx-1] = T_rp1

    # Neumann RP
    for i in range(Nx): #spodaj
        #A[i, i] = 1/(2*dy) #centralna diferenčna shema okoli 1.vrste - vzame 0. in 2. vrsto
        #A[i, i+2*Nx] = -1/(2*dy)
        A[i, i] = 1/(dy) #diferenčna shema naprej okoli 0.vrste - vzame 0. in 1. vrsto
        A[i, i+Nx] = -1/(dy)
        b[i] = 2 #+tok noter, -ven
    for j in range(1,Ny): #levo SYM RP - izolirano
        #A[j*Nx, j*Nx] = 1/(2*dx) #centralna diferenčna shema okoli 1.solpca - vzame 0. in 2. stolpec
        #A[j*Nx, j*Nx+2] = -1/(2*dx)
        A[j*Nx, j*Nx] = 1/(dx) #diferenčna shema naprej okoli 0.stolpca - vzame 0. in 1. stolpec
        A[j*Nx, j*Nx+1] = -1/(dx)
        b[j*Nx] = 0

    #desno zgoraj
    if insulated:
        for j in range(Ny//2,Ny): #desno (zgoraj) toplotno izolirano, NAMESTO T_rp
            #A[j*Nx+Nx-1, j*Nx+Nx-2] = 1/(2*dy) #centralna diferenčna shema okoli zadnjga stolpca - vzame -1. in -3. stolpec
            #A[j*Nx+Nx-1, j*Nx+Nx-1] = -1/(2*dy)
            A[j*Nx+Nx-1, j*Nx+Nx-2] = -1/(dy) #diferenčna shema nazaj okoli zadnjega stolpca - vzame -1. in -2. stolpec
            A[j*Nx+Nx-1, j*Nx+Nx-1] = 1/(dy)
            b[j*Nx+Nx-1] = 0
    else:
        for j in range(Ny//2,Ny):
            A[j*Nx+Nx-1, j*Nx+Nx-1] = 1
            b[j*Nx+Nx-1] = T_rp2 #R
    return A,b

def A_matrix(N, cylinder=True, insulated=False):
    Nx=Ny=N
    r=1 #polmer cevi - x
    L=1 #dolžina cevi - y
    dx=L/(Nx-1)
    dy=r/(Ny-1)

    # DE, RP
    A = np.zeros(((Nx)*(Ny), (Nx)*(Ny)))
    b = np.zeros(((Nx)*Ny))

    # DE
    for j in range(1,Ny-1):
        for i in range(1,Nx-1):
            ij_index = j*(Nx) + i #- Nx-1
            ij = 2*(1/dx**2 + 1/dy**2)
            ij_dx = -1/dx**2
            if cylinder: ij_dx += -1/(2*i*dx**2) #polarni dodatek!
            ij_dy = -1/dy**2

            #center
            A[ij_index, ij_index] += ij
            # okolica
            A[ij_index, ij_index-1] += ij_dx #levo
            A[ij_index, ij_index+1] += ij_dx #desno
            A[ij_index, ij_index-Nx] += ij_dy #dol
            A[ij_index, ij_index+Nx] += ij_dy #gor
            
            b[ij_index] = 0 #ni izvirov

    return RP(A,b,insulated,N)

#main solvers
def jacobi(A,q,N=1, insulated=False):
    Nx=Ny=int(len(q)**0.5)
    dx=dy=1./(Nx-1)

    # matrično
    # D_matrix = np.diagflat(D)
    # B = np.linalg.inv(D_matrix)
    # R = A - D_matrix

    D = np.diag(A)
    u = np.zeros_like(q)
    u_n = np.zeros_like(u)

    #iteracije
    for n in range(N):
        u_n = np.zeros_like(u)
        # Dirichlet RP
        for i in range(Nx): #gor
            u_n[-Nx+i] = q[-Nx+i]
        for j in range(Ny//2): #desno dol
            u_n[j*Nx+Nx-1] = q[j*Nx+Nx-1]
        for j in range(Ny//2,Ny-1): #desno gor
            u_n[j*Nx+Nx-1] = q[j*Nx+Nx-1]
        # Neumann RP
        for i in range(Nx): #spodaj
            u_n[i] = q[i]*dy + u[i+Nx]
        for j in range(1,Ny): #levo
            u_n[j*Nx] = q[j*Nx]*dx + u[j*Nx+1]
        if insulated:
            for j in range(Ny//2,Ny): #desno
                u_n[j*Nx+Nx-1] = q[j*Nx+Nx-1]*dx + u[j*Nx+Nx-2]
        
        for j in range(1,Ny-1):
            for i in range(1,Nx-1):
                ij_index = j*Nx + i
                u_n[ij_index] = (q[ij_index] + 1/dx**2*(u[ij_index+1] + u[ij_index-1]) +
                                               1/dy**2*(u[ij_index+Nx] + u[ij_index-Nx])) / D[ij_index]
        #u_n = np.dot(B, q-R@u) #matrično, hitreje za male Nx,Ny
        u = u_n
    return u
def gauss_seidel(A,q,N=1,insulated=False):
    Nx=Ny=int(len(q)**0.5)
    dx=dy=1./(Nx-1)

    D = np.diag(A)
    u = np.zeros_like(q)
    u_n = np.zeros_like(u)

    #iteracije
    for n in range(N):
        u_n = np.zeros_like(u)
        # Dirichlet RP
        for i in range(Nx): #gor
            u_n[-Nx+i] = q[-Nx+i]
        for j in range(Ny//2): #desno dol
            u_n[j*Nx+Nx-1] = q[j*Nx+Nx-1]
        for j in range(Ny//2,Ny-1): #desno gor
            u_n[j*Nx+Nx-1] = q[j*Nx+Nx-1]
        # Neumann RP
        for i in range(Nx): #spodaj
            u_n[i] = q[i]*dy + u[i+Nx]
        for j in range(1,Ny): #levo
            u_n[j*Nx] = q[j*Nx]*dx + u[j*Nx+1]
        if insulated:
            for j in range(Ny//2,Ny): #desno
                u_n[j*Nx+Nx-1] = q[j*Nx+Nx-1]*dx + u[j*Nx+Nx-2]
        for j in range(1,Ny-1):
            for i in range(1,Nx-1):
                ij_index = j*Nx + i
                u_n[ij_index] = 1/D[ij_index] * (q[ij_index] + 1/dx**2*(u[ij_index+1] + u_n[ij_index-1]) 
                                                             + 1/dy**2*(u[ij_index+Nx] + u_n[ij_index-Nx]))
        u = u_n
        #print(u)
    return u
def SOR(A,q,omega,N=1,insulated=False):
    Nx=Ny=int(len(q)**0.5)
    dx=dy=1./(Nx-1)

    D = np.diag(A)
    u_n = np.zeros_like(q)
    u = np.copy(u_n)

    #iteracije
    for n in range(N):
        #print(f'{n+1}/{N}',end='\r')
        # Dirichlet RP
        for i in range(Nx): #gor
            u_n[-Nx+i] = q[-Nx+i]
        for j in range(Ny//2): #desno dol
            u_n[j*Nx+Nx-1] = q[j*Nx+Nx-1]
        for j in range(Ny//2,Ny-1): #desno gor
            u_n[j*Nx+Nx-1] = q[j*Nx+Nx-1]
        # Neumann RP
        for i in range(Nx): #spodaj
            u_n[i] = q[i]*dy + u[i+Nx]
        for j in range(1,Ny): #levo
            u_n[j*Nx] = q[j*Nx]*dx + u[j*Nx+1]
        if insulated:
            for j in range(Ny//2,Ny): #desno
                u_n[j*Nx+Nx-1] = q[j*Nx+Nx-1]*dx + u[j*Nx+Nx-2]

        # u_n = SOR_i(D,q,omega,u_n)
        for j in range(1,Ny-1):
            for i in range(1,Nx-1):
                ij_index = j*Nx + i
                #u_gs = GS_i(A,q,u)
                u_gs = 1/D[ij_index] * (q[ij_index] + 1/dx**2*(u_n[ij_index+1] + u_n[ij_index-1]) 
                                                    + 1/dy**2*(u_n[ij_index+Nx] + u_n[ij_index-Nx]))
                u_n[ij_index] = u_n[ij_index] + omega*(u_gs - u_n[ij_index])

        u=np.copy(u_n)
    return u

#one iteration solver
def jac_i(D,q,u):
    Nx=Ny=int(len(q)**0.5)
    dx=dy=1./(Nx-1)

    u_n = np.copy(u)
    # Dirichlet RP
    for i in range(Nx): #gor
        u_n[-Nx+i] = q[-Nx+i]
    for j in range(Ny//2): #desno dol
        u_n[j*Nx+Nx-1] = q[j*Nx+Nx-1]
    for j in range(Ny//2,Ny-1): #desno gor
        u_n[j*Nx+Nx-1] = q[j*Nx+Nx-1]
    # Neumann RP
    for i in range(Nx): #spodaj
        u_n[i] = q[i]*dy + u[i+Nx]
    for j in range(1,Ny): #levo
        u_n[j*Nx] = q[j*Nx]*dx + u[j*Nx+1]

    for j in range(1,Ny-1):
        for i in range(1,Nx-1):
            ij_index = j*Nx + i
            u_n[ij_index] = (q[ij_index] + 1/dx**2*(u[ij_index+1] + u[ij_index-1]) + 
                           1/dy**2*(u[ij_index+Nx] + u[ij_index-Nx])) / D[ij_index]
    return u_n
def GS_i(D,q,u):
    Nx=Ny=int(len(q)**0.5)
    dx=dy=1./(Nx-1)

    u_n = np.copy(u)
    # Dirichlet RP
    for i in range(Nx): #gor
        u_n[-Nx+i] = q[-Nx+i]
    for j in range(Ny//2): #desno dol
        u_n[j*Nx+Nx-1] = q[j*Nx+Nx-1]
    for j in range(Ny//2,Ny-1): #desno gor
        u_n[j*Nx+Nx-1] = q[j*Nx+Nx-1]
    # Neumann RP
    for i in range(Nx): #spodaj
        u_n[i] = q[i]*dy + u[i+Nx]
    for j in range(1,Ny): #levo
        u_n[j*Nx] = q[j*Nx]*dx + u[j*Nx+1]

    for j in range(1,Ny-1):
        for i in range(1,Nx-1):
            ij_index = j*Nx + i
            u_n[ij_index] = 1/D[ij_index] * (q[ij_index] + 1/dx**2*(u[ij_index+1] + u_n[ij_index-1]) 
                                                         + 1/dy**2*(u[ij_index+Nx] + u_n[ij_index-Nx]))
    return u_n
def SOR_i(D,q,omega,u):
    Nx=Ny=int(len(q)**0.5)
    dx=dy=1./(Nx-1)

    # Dirichlet RP
    for i in range(Nx): #gor
        u[-Nx+i] = q[-Nx+i]
    for j in range(Ny//2): #desno dol
        u[j*Nx+Nx-1] = q[j*Nx+Nx-1]
    for j in range(Ny//2,Ny-1): #desno gor
        u[j*Nx+Nx-1] = q[j*Nx+Nx-1]
    # Neumann RP
    for i in range(Nx): #spodaj
        u[i] = q[i]*dy + u[i+Nx]
    for j in range(1,Ny): #levo
        u[j*Nx] = q[j*Nx]*dx + u[j*Nx+1]

    for j in range(1,Ny-1):
        for i in range(1,Nx-1):
            ij_index = j*Nx + i
            #u_gs = GS_i(A,q,u)
            u_gs = 1/D[ij_index] * (q[ij_index] + 1/dx**2*(u[ij_index+1] + u[ij_index-1]) 
                                                + 1/dy**2*(u[ij_index+Nx] + u[ij_index-Nx]))
            u[ij_index] = u[ij_index] + omega*(u_gs - u[ij_index])
    return u


def SOR_i2(D,q,omega,u_n):
    Nx=Ny=int(len(q)**0.5)
    dx=dy=1./(Nx-1)
    u = np.copy(u_n)

    # Dirichlet RP
    for i in range(Nx): #gor
        u_n[-Nx+i] = q[-Nx+i]
    for j in range(Ny//2): #desno dol
        u_n[j*Nx+Nx-1] = q[j*Nx+Nx-1]
    for j in range(Ny//2,Ny-1): #desno gor
        u_n[j*Nx+Nx-1] = q[j*Nx+Nx-1]
    # Neumann RP
    for i in range(Nx): #spodaj
        u_n[i] = q[i]*dy + u[i+Nx]
    for j in range(1,Ny): #levo
        u_n[j*Nx] = q[j*Nx]*dx + u[j*Nx+1]

    # u_n = SOR_i(D,q,omega,u_n)
    for j in range(1,Ny-1):
        for i in range(1,Nx-1):
            ij_index = j*Nx + i
            #u_gs = GS_i(A,q,u)
            u_gs = 1/D[ij_index] * (q[ij_index] + 1/dx**2*(u_n[ij_index+1] + u_n[ij_index-1]) 
                                                + 1/dy**2*(u_n[ij_index+Nx] + u_n[ij_index-Nx]))
            u_n[ij_index] = u_n[ij_index] + omega*(u_gs - u_n[ij_index])
    return u_n