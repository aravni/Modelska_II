#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
8.naloga: Galerkinov nastavek za reševanje problema lastnih vrednosti na polkrožni domeni.
Povzeto po github [vir]
"""

import numpy as np
import scipy.special as sp
import scipy.linalg as la

import matplotlib.pyplot as plt
from matplotlib import rc

from timeit import default_timer as timer


# Plotting setup
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

import matplotlib
matplotlib.use('TkAgg')

# trial functions
def Psi(m, n, xi, phi, ftype='maii'):
	if ftype == 'mfp':
		return xi**(2*m+1) * (1-xi)**(n+1) * np.sin((2*m+1)*phi)
	elif ftype == 'maii':
		return xi**(m+n) * (1-xi) * np.sin(m*phi)

# matrix A
def A(nt, m, n, ftype='maii'):
	if ftype == 'mfp':
		return -np.pi/2 * ((n+1)*(nt+1)*(4*m+3))/(n+nt+4*m+4) * sp.beta(n+nt+1, 4*m+3)
	elif ftype == 'maii':
		return np.pi/2 * ((n*nt)/(2*m+n+nt)-(2*n*nt+n+nt)/(2*m+n+nt+1)+((n+1)*(nt+1))/(2*m+2+n+nt))
		#return np.pi/2 * ((m+n)*(1-m-n)/(2*m-1+n+nt)+2*(m+n)**2/(2*m+n+nt)-(m+n+2*m*n+n**2)/(2*m+1+n+nt)-2*m**2/(2*m+2+n+nt)+m**2/(2*m+3+n+nt))

# matrix B
def B(nt, m, n, ftype='maii'):
	if ftype == 'mfp':
		return 0 # WIP
	elif ftype == 'maii':
		return np.pi/2 * (1/(2*m+2+n+nt)-2/(2*m+3+n+nt)+1/(2*m+4+n+nt))
		#return np.pi/2 * (2/((2*m+1+n+nt)*(2*m+2+n+nt)*(2*m+3+n+nt)))

# do the calculation
def calc(M, N, eigs, to_calc_u=True, to_plot_u=False, xi_num=1000, phi_num=1000, Nc=2,Nr=2):
	t1=timer()
	# Setup of parameters for calculating u
	if to_calc_u:
		
		# Mesh grid for plotting and calculating u
		# ksi is wighted more towards outer edge
		xi = np.linspace(0, 1, xi_num)
		phi = np.linspace(0, 1, phi_num)**(1/3) * np.pi

		r, theta = np.meshgrid(xi, phi)

		# Initial matrix of values and error
		u = np.zeros((M, eigs, xi_num, phi_num))
		deltaumax_lst = []
		lmd = np.zeros((M, eigs))

	# size of vector set
	m_lst = range(1, M)
	n_lst = range(N)

	# Am=np.zeros((N*M, N*M))
	# Bm=np.zeros((N*M, N*M))
	# system is decomposable in m, so we solve for every m separately
	for m in m_lst:

		# indicator
		# print(f'Solving for m={m} out of M={M} with N={N}')

		# for each m we solve an eigenvalue problem
		# with matrix A
		A_mat = np.fromfunction(lambda nt, n: A(nt, m, n), (N, N))

		# and matrix B
		B_mat = np.fromfunction(lambda nt, n: B(nt, m, n), (N, N))
		
		# Am[(m-1)*N:(m)*N, (m-1)*N:(m)*N] = A_mat
		# Bm[(m-1)*N:(m)*N, (m-1)*N:(m)*N] = B_mat

		# solve the system
		ep, vp = la.eigh(A_mat, B_mat, subset_by_index=[0, eigs-1])
		# print(np.sqrt(ep))
	# plt.matshow(Am)
	# plt.colorbar()
	# plt.show()

	# # solve the system
	# ep, vp = la.eig(Am, Bm)

	# for m in m_lst:
		if to_calc_u:
		
			# vector of matrices for trial functions
			psi = np.fromfunction(lambda n: Psi(m, n, xi[:, None, None],
				phi[:, None]), (N,))

			# print(psi.shape, vp.shape)

			# build u
			for i in range(eigs):
				
				# part of solution up
				up = np.sum(psi * vp[:,i].reshape((1,1,-1)), axis=2)
				deltaumax_lst.append(np.max(np.abs(up)))
				u[m,i] += up
		
			lmd[m,:] = ep
	lmd = np.concatenate(lmd)
	k = sorted(np.sqrt(lmd[lmd > 0]))
	u = u.transpose(1, 0, 3,2)
	# print(k)
	t2=timer()

	# Plotting velocity profile
	if to_plot_u and to_calc_u:
		# n_eig = 2
		# fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
		# ax.set_title(f"$k_{{{n_eig}}}={k[n_eig-1]:.5g}$")
		# ax.set_xlim(0, np.pi)
		# # plt.contour(theta, r, u[2,0].T, 7, colors='k',
		# # 	linestyles='dashed', linewidths=1)
		# surf = ax.contourf(theta, r, u[0,n_eig], 10, cmap='turbo')
		# clb = fig.colorbar(surf, label='$u$')
		# plt.show()

		X=r*np.cos(theta)
		Y=r*np.sin(theta)
		vec=u[0]
		fig = plt.figure(figsize=(2*Nc*2+0.5, Nr*2))
		gs = fig.add_gridspec(Nr,Nc)
		ax = gs.subplots(sharex=True, sharey=True)
		vmax,vmin=np.max(vec[:Nr*Nc]), np.min(vec[:Nr*Nc])
		i=j=0
		m,n=1,1
		for l in range(Nr*Nc):
			cs=ax[j,i].contourf(X,Y, vec[l], vmin=vmin, vmax=vmax, levels=10, cmap='turbo')
			ax[j,i].grid(False)
			ax[j,i].axis('off')
			ax[j,i].set_xticks([])
			ax[j,i].set_yticks([])
			ax[j,i].set_title(f'$k_{{{i+1}}}$ = {k[l]:.3g}', fontsize=11)
			i+=1
			n+=1
			if n+1>m: n=1; m+=1
			if i==Nc: i=0; j+=1
		cbar_ax = fig.add_axes([0.86, 0, 0.07, 1])
		cbar_ax.axis('off')
		cbar=fig.colorbar(cs, ax=cbar_ax, orientation='vertical')
		cbar.ax.tick_params(labelsize=11)
		cbar.set_label('u', fontsize=12)
		plt.show()

	if to_calc_u:
		return k, u, t2-t1
	else: return None

#------------------------------------------------------------------------------
# Plotting the profile
if __name__ == '__main__':
	# Parameters M and N
	M = 7
	N = 10

	# calculate C and plot u 
	lmd, _, _ = calc(M, N, 5, to_calc_u=True, to_plot_u=True)
#------------------------------------------------------------------------------
 
