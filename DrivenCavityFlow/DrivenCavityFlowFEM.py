"""
Program computes the driven cavity flow by using the finite element
method in fenics. This is a standalone script made as simple as possible
for educational purposes
"""
from fenics import *

# Set up the space
mesh = RectangleMesh(Point(0,0), Point(1,1), 3, 3)
VV = VectorFunctionSpace(mesh, "CG", 1)
V = FunctionSpace(mesh, "CG", 1)
uvv = TrialFunction(VV)
vvv = TestFunction(VV)
uv = TrialFunction(V)
vv = TestFunction(V)

# Set vectors
up = Function(VV)  # Old velocity vector
un = Function(VV)  # New velocity vector
rp = Function(V)   # Old density
rn = Function(V)   # New density

# Initial condition
up.interpolate(Constant((0,0)))
rp.interpolate(Constant(0))

# Boundary functions
def lid(x, on): return on and abs(x[1] - 1) < DOLFIN_EPS
def left(x, on): return on and abs(x[0]) < DOLFIN_EPS
def right(x, on): return on and abs(x[0] - 1) < DOLFIN_EPS
def bottom(x, on): return on and abs(x[1]) < DOLFIN_EPS

# Create Boundary conditions
lidbc = DirichletBC(VV, Constant((1,0)), lid)
leftbc = DirichletBC(VV, Constant((0,0)), left)
rightbc = DirichletBC(VV, Constant((0,0)), right)
bottombc = DirichletBC(VV, Constant((0,0)), bottom)

bcs = [lidbc, leftbc, rightbc, bottombc]

# Set boundary functions to the initial condition
[bc.apply(up.vector()) for bc in bcs]
print up.vector().array()

# Begin solver
