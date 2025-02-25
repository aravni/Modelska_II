#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 6. naloga: PDE lastne rešitve
# homogena opna

import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer
from tqdm import tqdm
from numba import njit, jit
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from scipy import sparse
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor

plt.style.use('aljaz')
np.set_printoptions(linewidth=200, precision=3, suppress=True)

def is_diag(A):
    """Check if A is near-diagonal by comparing off-diagonal to diagonal elements.
    If the result is small, it's near-diagonal."""
    A_off_diag = A.copy()
    np.fill_diagonal(A_off_diag, 0)  # Zero out the diagonal elements
    off_diag_norm = np.linalg.norm(A_off_diag, ord='fro')  # Frobenius norm of off-diagonal elements
    diag_norm = np.linalg.norm(np.diag(A), ord=2)  # Norm of diagonal elements
    return off_diag_norm / diag_norm if diag_norm > 0 else np.inf
@njit(cache=True)
def max_el(A):
    '''največji izven-diagonalni element simetrične matrike'''
    N=len(A)
    A_max=0.
    I=J=0
    for i in range(N):
        for j in range(i+1,N): #sym
            if np.abs(A[i,j]) > A_max:
                A_max = np.abs(A[i,j])
                I,J = i,j
    return A_max,I,J
@njit(cache=True)
def jacobi(A, tol=1e-6):
    '''Jacobi metoda za diagonalizacijo simetrične matrike A.
    Vrne diagonalne elemente matrike A (lastne vrednosti) in matriko J (lastni vektorji).'''
    N=A.shape[0]
    Nmax=5*N**2 # da skonvergira
    J = np.eye(N) #Jakobijeve rotacije
    for n in (range(Nmax)):
        A_max,j,k = max_el(A)
        if A_max < tol: break
        if A[j, j] == A[k, k]: theta = np.pi/4 #singularno
        else: theta = 0.5*np.arctan(2*A[j,k]/(A[j,j]-A[k,k]))
        c=np.cos(theta)
        s=np.sin(theta)
        # R = np.array([[c,-s],[s, c]]) #rotacijska matrika
        # ind=np.ix_([j,k],[j,k])
        # J[ind] += R
        # A[ind] = R.T @ A[ind] @ R #iteracija NAPAKA?!
        # A = J.T @ A @ J

        # A[:, (j, k)] = A[:, (j, k)] @ R
        # A[(j, k), :] = R.T @ A[(j, k), :]
        # J[:, (j, k)] = J[:, (j, k)] @ R
        # numba njit ne podpira python indexiranja
        for i in range(N):
            A_ij, A_ik = A[i, j], A[i, k]
            A[i, j] = c * A_ij + s * A_ik
            A[i, k] = -s * A_ij + c * A_ik
        for i in range(N):
            A_ji, A_ki = A[j, i], A[k, i]
            A[j, i] = c * A_ji + s * A_ki
            A[k, i] = -s * A_ji + c * A_ki
        for i in range(N):
            J_ij, J_ik = J[i, j], J[i, k]
            J[i, j] = c * J_ij + s * J_ik
            J[i, k] = -s * J_ij + c * J_ik
    if A_max > tol: print('Jacobi did not converge in', Nmax, 'iterations:', A_max, '>', tol) #@jit brez f stringa
    else: print(f'Jacobi converged in', n, 'iterations:', A_max, '<', tol)
    # if A_max > tol: print(f'Jacobi did not converge in {Nmax} iterations: {A_max:.1E} > {tol:.1E}!')
    # else: print(f'Jacobi converged in {n} iterations: {A_max:.1E} < {tol:.1E}.')
    # print('test diagonalnosti: ', round(is_diag(A),2)) #ni optimizirano!
    ind=np.argsort(np.diag(A))
    w=np.diag(A)[ind] #lastne vrednosti
    x=J[:,ind] #lastni vektorji
    return w,x,A

# matrika A
def mat(N,h):
    A = np.zeros((N**2,N**2))
    for i in range(N):
        for j in range(N):
            # if i==0 or i==N-1 or j==0 or j==N-1: #RP?
            #     A[i+j*N,i+j*N] = -1#-4/h/h
            #     pass
            # else:
            #     A[i+j*N,i+j*N] = -4/h/h
            #     A[i+j*N,i+(j-1)*N] = 1/h/h
            #     A[i+j*N,i+(j+1)*N] = 1/h/h
            #     A[i+j*N,i-1+j*N] = 1/h/h
            #     A[i+j*N,i+1+j*N] = 1/h/h
            A[i+j*N, i+j*N] = 4/h/h  # Main diagonal
            if i > 0:
                A[i+j*N, i-1+j*N] = -1/h/h  # Lower diagonal
            if j > 0:
                A[i+j*N, i+(j-1)*N] = -1/h/h  # Left diagonal
            if i < N-1:
                A[i+j*N, i+1+j*N] = -1/h/h  # Upper diagonal
            if j < N-1:
                A[i+j*N, i+(j+1)*N] = -1/h/h  # Right diagonal
    return A

# problem last vrednosti
def last(A):
    # w1, x1 = np.linalg.eig(A)
    t1=timer()
    w, x = np.linalg.eigh(A) #sym -> hitreje!
    t_eig=timer()-t1
    print('t eig = ', round(t_eig,3), ' s')
    # print(np.sum(np.abs(w1-w)), np.sum(np.abs((x1-x)))) #ok za nizke lastne vrednosti
    return w, x, t_eig
def sp(A, Ns):
    t1=timer()
    # ws, xs = sparse.linalg.eigsh(A, k=Ns, which='LM') #velike vrednosti -> OK
    ws, xs = sparse.linalg.eigsh(A, k=Ns, which='SM') #male vrednosti -> slabše!?
    # ws, xs = sparse.linalg.eigsh(A, k=Ns, sigma=0, which='LM') #shift invert mode za male vrednosti
    t_sp=timer()-t1
    # print(ws.shape, xs.shape, w.shape, x.shape)
    print('t sparse = ', round(t_sp,3), ' s')
    return ws, xs, t_sp
def jac(A):
    t1=timer()
    wj, xj, Aj = jacobi(A, tol=1e-6)
    t_j=timer()-t1
    print('t jacobi = ', round(t_j,3), ' s')
    print('test diagonalnosti: ', is_diag(Aj))
    return wj, xj, t_j
# print(is_diag(np.eye(N**2)+1e-3*np.random.random((N*N,N*N))), 'test')
# print(is_diag(np.diag(np.random.random(N**2)*1e1)), 'test')
def parallel(N,a=1,N_sparse=9):
    h=a/(N-1)
    A=mat(N,h)
    w,x,t_eig=last(A)
    ws,xs,t_sp=sp(A, N_sparse)
    if N<50: wj,xj,t_j=jac(A)
    else: t_j=None #predolgo!
    return (t_eig,t_sp,t_j)
if __name__ == '__main__':
    # N_list=np.linspace(3,100,18,dtype=int) #test

    # # Function to update the progress bar
    # def update_progress(*args):
    #     progress_bar.update()
    # # Initialize the progress bar
    # progress_bar = tqdm(total=len(N_list), desc='Processing')
    # progress_bar.update()

    # with ProcessPoolExecutor(max_workers=max([len(N_list),18])) as executor:
    #     # Create a list of futures
    #     futures = [executor.submit(parallel, N) for N in N_list]
    #     for future in futures:
    #         future.add_done_callback(update_progress)
    #     # Wait for all futures to complete and collect results
    #     results = [future.result() for future in futures]
    #     # Close the progress bar
    #     progress_bar.close()
    # čas = {'eig':[], 'sparse':[], 'jacobi':[]}
    # for result in results:
    #     čas['eig'].append(result[0])
    #     čas['sparse'].append(result[1])
    #     čas['jacobi'].append(result[2])

    # plt.plot(N_list, čas['eig'], '.-', label='scipy.eigh')
    # plt.plot(N_list, čas['sparse'], '.-', label='ARPACK')
    # plt.plot(N_list, čas['jacobi'], '.-', label='Jacobi')
    # # plt.xscale('log')
    # plt.yscale('log')
    # plt.xlabel('$N$')
    # plt.ylabel('$t$ [s]')
    # plt.legend()
    # # plt.savefig('./slike/cas_homog.pdf')
    # plt.show()

    ##############################################################
    N=10 #št. računskih točk po osi
    a=b=1 #širina=višina
    h=a/(N-1)

    # snovne lastnosti
    gamma=1 # RP sila vpetja, N/m
    d=1 #debelina, m
    rho=1 #kg/m^3
    c=np.sqrt(gamma/rho/d) #hitrost valovanja, m/s
    print('c = ',c,' m/s')

    Nc,Nr=3,3 #št. stoplcev, vrstic
    l=2

    A = mat(N,h)
    w,x,t_eig = last(A)

    fig=plt.figure(figsize=(Nc*l,Nr*l+0.5))
    gs=fig.add_gridspec(Nr,Nc, hspace=0.11, wspace=0.11)
    ax=gs.subplots(sharex=True, sharey=True)
    vmax,vmin=np.max(x[:,:Nr*Nc]), np.min(x[:,:Nr*Nc])
    i=j=0
    for k in range(Nr*Nc):
        ax[j,i].imshow((x[:,k]).reshape(N,N), origin='lower', extent=(0,a,0,b), aspect='equal', vmax=vmax, vmin=vmin)
        ax[j,i].grid(False)
        ax[j,i].set_xticks([])
        ax[j,i].set_yticks([])
        # ax[j,i].set_title(f'$\lambda_{{{k+1}}}$', fontsize=9)
        ax[j,i].set_title(f'$f_{{{k+1}}} = {(np.sqrt(w[k])*c/2/np.pi):<.2f} Hz$', fontsize=9)
        i+=1
        if i==Nc: i=0; j+=1
    cbar_ax = fig.add_axes([0.86, 0, 0.07, 1])
    cbar_ax.axis('off')
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    sm = cm.ScalarMappable(cmap='turbo', norm=norm)
    sm.set_array([])
    cbar=fig.colorbar(sm, ax=cbar_ax)
    cbar.ax.tick_params(labelsize=11)
    # cbar.set_ticks(np.arange(vmin,vmax,0.01))
    # plt.tight_layout()
    # plt.savefig('./slike/lastne_homog.pdf')
    plt.show()

    # sparse
    ws,xs,t_sp = sp(A, Nc*Nr)
    fig=plt.figure(figsize=(Nc*l,Nr*l+0.5))
    gs=fig.add_gridspec(Nr,Nc, hspace=0.11, wspace=0.11)
    ax=gs.subplots(sharex=True, sharey=True)
    vmax,vmin=np.max(np.abs(xs[:,:Nr*Nc])), np.min(np.abs(xs[:,:Nr*Nc]))
    vmax,vmin=np.max((xs[:,:Nr*Nc])), np.min((xs[:,:Nr*Nc]))
    i=j=0
    for k in range(Nr*Nc):
        # ax[j,i].imshow(np.abs(xs[:,k]).reshape(N,N), origin='lower', extent=(0,a,0,b), aspect='equal', vmax=vmax, vmin=vmin)
        ax[j,i].imshow((xs[:,k]).reshape(N,N), origin='lower', extent=(0,a,0,b), aspect='equal', vmax=vmax, vmin=vmin)
        ax[j,i].grid(False)
        ax[j,i].set_xticks([])
        ax[j,i].set_yticks([])
        ax[j,i].set_title(f'$\lambda_{{{k+1}}}$', fontsize=9)
        # ax[j,i].set_title(f'$\lambda_{{{k+1}}} = {(np.sqrt(-ws[k])*c/2/np.pi):<.2f} Hz$', fontsize=9)
        i+=1
        if i==Nc: i=0; j+=1
    cbar_ax = fig.add_axes([0.86, 0, 0.07, 1])
    cbar_ax.axis('off')
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    sm = cm.ScalarMappable(cmap='turbo', norm=norm)
    sm.set_array([])
    cbar=fig.colorbar(sm, ax=cbar_ax)
    cbar.ax.tick_params(labelsize=11)
    # cbar.set_ticks(np.arange(vmin,vmax,0.01))
    # plt.tight_layout()
    # plt.savefig('./slike/sparse_homog.pdf')
    plt.show()

    # jacobi
    wj,xj,t_j = jac(A)
    fig=plt.figure(figsize=(Nc*l,Nr*l+0.5))
    gs=fig.add_gridspec(Nr,Nc, hspace=0.11, wspace=0.11)
    ax=gs.subplots(sharex=True, sharey=True)
    vmax,vmin=np.max((xj[:,:Nr*Nc])), np.min((xj[:,:Nr*Nc]))
    i=j=0
    for k in range(Nr*Nc):
        ax[j,i].imshow((xj[:,k]).reshape(N,N), origin='lower', extent=(0,a,0,b), aspect='equal', vmax=vmax, vmin=vmin)
        ax[j,i].grid(False)
        ax[j,i].set_xticks([])
        ax[j,i].set_yticks([])
        ax[j,i].set_title(f'$\lambda_{{{k+1}}}$', fontsize=9)
        # ax[j,i].set_title(f'$\lambda_{{{k+1}}} = {(np.sqrt(-wj[k])*c/2/np.pi):<.2f} Hz$', fontsize=9)
        i+=1
        if i==Nc: i=0; j+=1
    cbar_ax = fig.add_axes([0.86, 0, 0.07, 1])
    cbar_ax.axis('off')
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    sm = cm.ScalarMappable(cmap='turbo', norm=norm)
    sm.set_array([])
    cbar=fig.colorbar(sm, ax=cbar_ax)
    cbar.ax.tick_params(labelsize=11)
    # cbar.set_ticks(np.arange(vmin,vmax,0.01))
    # plt.tight_layout()
    # plt.savefig('./slike/jacobi_homog.pdf')
    plt.show()

    plt.plot(np.sqrt(w)*c/2/np.pi, 'o', label='scipy.eigh')
    plt.plot(np.sqrt(ws)*c/2/np.pi, 'x', label='ARPACK')
    plt.plot(np.sqrt(wj)*c/2/np.pi, '+', label='Jacobi')
    plt.legend()
    plt.xlabel('indeks')
    plt.ylabel('$f$ [Hz]')
    plt.minorticks_on()
    # plt.savefig('./slike/frekvence_homog.pdf')
    plt.show()