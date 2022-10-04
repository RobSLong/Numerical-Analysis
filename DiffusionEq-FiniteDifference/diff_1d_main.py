# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 22:12:26 2022

@author: LongR
"""

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.family'] = 'serif'

def diffusion_ftcs(T0, nt, dt, dx, alpha, BC_type = ['Dirichlet', 'Dirichlet']):
    """
    Function to numerically solve diffusion equation using 
    FTCS Finite Difference scheme
    
    Parameters
    ----------
    T0 : Initial condition for temperature (1D array of floats)
    nt : Number of time steps (integer)
    dt : Time step size (integer)
    dx : Spatial grid spacing (integer)
    alpha : Thermal diffusivity of rod, m^2/s (float)
    BC_type : Left and right boundary conditions (string) The default is ['Dirichlet', 'Dirichlet'].

    Returns
    -------
    T : Temperature along the domain (1D array of floats)
    Tplot : Temperature along the domain saved sparsely for plotting (array of floats)
    """
    Tplot = []; # initialise array
    T = T0.copy() # 
    sigma = alpha * dt/dx**2  # coefficient in Finite Difference scheme
    for i in range(nt): # loop over time
        T[1:-1] = T[1:-1] + sigma * (T[2:] - 2*T[1:-1] + T[:-2])
        if BC_type[0] == 'Neumann': # implement left Neumann BC
            T[0] = T[1]
        if BC_type[1] == 'Neumann': # implement right Neumann BC
            T[-1] = T[-2]
        if i % 200 == 0: # save T data systematically for plotting
            Tplot.append(T.copy())
    return T, Tplot

def plotter(file_name):
    """
    Function to create and save solution as a png
    
    Parameter
    ---------
    file_name : name to save the output file (string)
    """
    # figure formatting
    plt.figure(figsize=(4,3.3), dpi=200)
    plt.grid(alpha=0.4)
    plt.xlabel('Horizontal position, $x$'); plt.ylabel('Temperature, $T$')
    plt.xlim([-0.02,1.02]); plt.ylim([-0.02,1.02])
    plt.plot(x,IC,'k-',lw=3)
    # create orange to purple colour scheme
    R = np.linspace(1,0.4,len(Tplot))
    B = np.linspace(0,1,len(Tplot))
    G = np.linspace(0.647,0,len(Tplot))
    plt.plot(x, np.exp(-1* np.pi **2 * alpha * (600*dt))*np.sin(np.pi*x),'ko')    
    for k in range(0,len(Tplot)):
        plt.plot(x,Tplot[k], color=[R[k], G[k], B[k]])
    
    plt.tight_layout()
    plt.savefig(file_name+ '.png')
    plt.show()

from scipy import linalg


def  lhs(N, sigma):
    D = np.diag((2 + 1/sigma) * np.ones(N))
    D[-1,-1] = 1 + 1/sigma
    u = np.diag(-1 * np.ones(N-1), k=1)
    L = np.diag(-1 * np.ones(N-1), k=-1)
    A = D + u + L
    return A

def rhs(T, sigma, qdx):
    b = T[1:-1] / sigma
    b[0] += T[0]
    b[-1] += qdx
    return b
    
def btcs_impl(T0, nt, dt, dx, alpha, q):
    sigma = alpha * dt/ dx**2
    A = lhs(len(T0) - 2, sigma)
    T = T0.copy()
    for i in range(nt):
        # Create right hand side of system
        b = rhs(T, sigma, q * dx)
        T[1:-1] = linalg.solve(A, b)
        # Apply RHS Neumann
        T[-1] = T[-2] + q * dx
    return T

#%% Case study 1
L = 1; nx = 51; 
dx = L / (nx - 1); x = np.linspace(0,L,nx)
alpha = 1.22e-3
q= 0.
sigma = 0.5
dt = sigma * dx**2 / alpha
nt = 100
#T0 = np.zeros(nx)
#T0[0] = 1.0
T0 = np.sin(np.pi*x)
T = btcs_impl(T0,nt,dt,dx,alpha, q)

plt.figure(figsize=(4,3.3), dpi=200)
plt.grid(alpha=0.4)
plt.xlabel('Horizontal position, $x$'); plt.ylabel('Temperature, $T$')
#plt.xlim([-0.02,1.02]); plt.ylim([-0.02,1.02])
plt.plot(x,T,'k-',lw=3)
plt.tight_layout()
plt.show()

#%% Case study 1
L = 1; nx = 51; 
dx = L / (nx - 1); x = np.linspace(0,L,nx)
alpha = 1.22e-3
nt = 1801


sigma = 0.5
dt = sigma * dx**2 / alpha

#%%

IC = np.sin(np.pi*x)
T, Tplot = diffusion_ftcs(IC, nt, dt, dx, alpha)

plotter('diff_sinIC')

#%% Case study 1
nt = 1201
IC = np.zeros(nx)
IC[0] = 1; IC[-1] = 0.2
T, Tplot = diffusion_ftcs(IC, nt, dt, dx, alpha)
plotter('diff_DirDir')
#%%
nt = 20001
IC = np.zeros(nx)
IC[0] = 1
IC[-1] = 0.2
T, Tplot = diffusion_ftcs(IC, nt, dt, dx, alpha, ['Dirichlet', 'Neumann'])
plotter('diff_DirNeu')