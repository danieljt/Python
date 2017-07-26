"""
Test program for running a driven cavity flow problem from the
DrivenCavityFlowFDM.py program. 
"""
from DrivenCavityFlowFDM import *

# Set constants and run solver
nx = 40
ny = 40
L = 1
H = 1
Re = 1
M = 0.1
dt = 0.0001
T = 0.1
U = 1
solver(nx,ny,H,L,dt,T,Re,M,U)
