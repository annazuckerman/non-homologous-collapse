import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import ScalarFormatter
import os
import argparse

# parse input arguments
parser = argparse.ArgumentParser(description='')
parser.add_argument('--epsilon',
                    dest='epsilon',
                    type=float,
                    help='Degree of initial overdensity.',
                    required=True)
parser.add_argument('--outdir',
                    dest='outdir',
                    type=str,
                    help='Directory to save output plots.',
                    default=True)
args = parser.parse_args()
epsilon = args.epsilon
outdir = args.outdir

def diff(x, y):
    '''
    Approximate the derivative of y with respect to x, using three-point Lagrangian interpolation.
    Following the method of the IDL DERIV method described here: www.nv5geospatialsoftware.com/docs/DERIV.html
    Inputs:
        x (numpy array): x values for differentiation
    Returns:
        y (numpy array): y values corresponding to x values
    '''
    dydx = np.ones(len(y)) * np.nan
    y0 = y[0]
    y1 = y[1]
    y2 = y[2]
    x01 = x[0] - x[1]
    x02 = x[1] - x[2]
    x12 = x[1] - x[2]
    dydx[0] = y0*(x01 + x02)/(x01*x02) - y1*x02/(x01*x12) + y2*x01/(x02*x12)
    for i in range(1,len(y)-1):
        y0 = y[i-1]
        y1 = y[i]
        y2 = y[i+1]
        x01 = x[i-1] - x[i]
        x02 = x[i-1] - x[i+1]
        x12 = x[i] - x[i+1]
        dydx[i] = y0*x12/(x01*x02) + y1*(1/x12 - 1/x01) - y2*x01/(x02*x12)
    y0 = y[-3]
    y1 = y[-2]
    y2 = y[-1]
    x01 = x[-3] - x[-2]
    x02 = x[-3] - x[-1]
    x12 = x[-2] - x[-1]
    dydx[-1] = -y0*x12/(x01*x02) + y1*x02/(x01*x12) - y2*(x02 + x12)/(x02*x12) 
    return dydx


# define constants
Gconst = 6.6732e-08 # universal gravitational constant
xMsun  = 1.989e+33 # mass of sun in g
xRsun  = 6.96e+10 # radius of sun in cm

# set up initial properties of the model
bigR0   = xRsun # initial stellar radius
rhoc0   = xMsun/(4*np.pi*(bigR0**3)/3) # initial central density

# define min and and max colapse times
tDmin = np.sqrt(3*np.pi/32/Gconst/rhoc0)
tDmax = tDmin / np.sqrt(1 - 0.75*epsilon)

# set up grid in r0, and define mass inside these initial shells
nr0  = 800 # number of shells
amin = 1e-3*bigR0 # 
amax = bigR0 - amin
r0 = np.linspace(amin,amax,nr0) 
Mr0  = 4*np.pi/3 * rhoc0 * (r0**3) * (1 - 0.75*epsilon*r0/bigR0) 

# set up "time" grid (initially in xi)
nt   = 900
xi   = np.linspace(0,nt-1,nt)/(nt-1)*(0.5*np.pi * 0.99999) 

# determine r(r0,xi) and t(r0,xi) in 2D grids
r = np.zeros([nr0,nt])
t = np.zeros([nr0,nt])
for j in range(nt):
    for i in range(nr0):
        r[i,j] = r0[i] * np.cos(xi[j])**2 
        t[i,j] = np.sqrt(r0[i]**3/(2*Gconst*Mr0[i])) * (xi[j] + 0.5*np.sin(2*xi[j]))        
tnorm = t/tDmax
rnorm = r/bigR0

# compute slices of constant time (new grid in t-prime)
# and numerically differentiate
nslice = 100
tslice = np.linspace(0,0.95*tDmin,nslice)
rslice = np.zeros([nslice,nr0])
for i in range(nr0):  
    rslice[:,i] = np.interp(tslice, t[i,:], r[i,:]) 
rho = np.zeros([nslice,nr0])
for j in range(nslice):
    x = rslice[j,:]
    y = Mr0
    dMdr = diff(x, y)
    rho[j,:] = dMdr / (4*np.pi*x[:]**2)   

# make plots
fs1 = 14
fs2 = 16
plt.figure(figsize = [30,13], dpi = 200)
plt.subplot(221)
for i in range(nr0):
    if (i%100 == 0) | (i == 799): # for clarity in plotting
        l1, = plt.plot(tnorm[i,:], rnorm[i,:], label = 'r0='+str(round((r0[i]/xRsun),3))+'R0') # r/R0 vs time/t_max 
        color = l1.get_color()
        plt.plot(tslice/tDmax, rslice[:,i]/bigR0, '.', color = color)
plt.plot([],[],'.', color = 'k', label = 'Interpolation points')
plt.xlabel(r't/t$_{max}$', fontsize = fs2)
plt.ylabel(r'r/R$_0$', fontsize = fs2)
plt.legend(fontsize = fs1)
plt.subplot(223)
for i in range(nr0):
    if (i%100 == 0) | (i== 799): # for clarity in plotting
        plt.plot(tnorm[i,:], r[i,:]/r0[i], label = 'r0='+str(round((r0[i]/xRsun),3))+'R0')
plt.xlabel(r't/t$_{max}$', fontsize = fs2)
plt.ylabel(r'r/r$_0$', fontsize = fs2)
plt.legend(fontsize = fs1)
plt.subplot(222)
for j in range(nslice):
    if (j%15 == 0) | (j== nslice - 1): # for clarity in plotting
        xpoly = rslice[j,:]/rslice[j,nr0-1] # rslice over final r0 rslice (divide just for normalizing)
        ypoly = rho[j,:]/rho[j,nr0-1] # rslice over final r0 rslice (divide just for normalizing)
        idxs = xpoly >= 0.05
        plt.plot(xpoly[idxs],ypoly[idxs],  label = r't/t$_{min}=$'+str(round(tslice[j]/tDmin,2)))
if epsilon == 0:
    plt.ylim([0.97,1.03])
else:
    plt.semilogy()
    plt.gca().yaxis.set_major_formatter(ScalarFormatter())
plt.xlabel(r'r/r$_0$', fontsize = fs2)
plt.ylabel(r'$\rho$ (normalized)', fontsize = fs2)
plt.legend(fontsize = fs1)
plt.subplot(224)
rho_c =  rho[:,10] 
rho_avg = xMsun / (4*np.pi*rslice[:,-1]**3/3) 
plt.plot(tslice/tDmax, rho_c/rho_avg)
plt.xlabel(r't/t$_{max}$', fontsize = fs2)
plt.ylabel(r'$\rho_c/\rho_{avg}$', fontsize = fs2)
if epsilon == 0:
    plt.ylim([0.97,1.03])
plt.tight_layout()
if not os.path.exists(outdir):
    os.mkdir(outdir)
    
outfile = 'non_homologous_plots_eps=' + str(epsilon) + '.png'
plt.savefig(outdir + '/' + outfile)
