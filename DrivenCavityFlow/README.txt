This folder holds programs for solving the lid driven driven cavity flow in a
rectangular domain. This is a well known and researched field of fluid
mechanics, and these programs are in no way unique when compared to the vast
amount of scripts that are already produced. Yet, this folder contains two
programs solving the same problem. DrivenCavityFlowFDM.py solves the cavity
flow using the explicit MacCormack scheme, while the DrivenCavityFlowFEM.py
solves the problem using the finite element method. Both programs use
non-dimentionalized values. For usage, see the MainFDM.py and MainFEM.py
for running the FDM and FEM scheme respectively. The theory behind the methods
is assumed to be familiar to the reader, but can be found at the following
sources:
Pijush K. Kundu, Ira M. Cohen, David R. Dowling Fluid mechanics 5th edition
www.fenicsproject.org

PACAKAGE VERSIONS:
python     2.7.12
numpy      1.11.0
matplotlib 1.5.1
fenics     2016.1.0
scitools


