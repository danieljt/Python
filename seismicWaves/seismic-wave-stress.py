from dolfin import *
import mayavi as ma
import math as mt

def solver(L,h,xel,yel,dt,T,lamda,mu,rho,A,omega,nx,ny,wavetype,viz,savefile):
    """
    Function for solving the seismic wave equation on a rectangular
    domain with with inhomogeneous dirichlet bcs on 3 sides, and with
    a given stress on the top. The solver
    is tested with either a p wave or s-wave solution.
    INPUT:
    L       : Length of domain
    h       : Height of domain
    xel     : Number of elements per unit length in x direction
    yel     : Number of elements per unit length in y direction
    dt      : Time step
    T       : Total simulation time
    lamda   : lamees first parameter
    mu      : shear modulus
    rho     : density of material
    A       : Amplitude of test solution
    omega   : angular frequency of test solution
    nx      : component of normal vector of test solution in x direction
    ny      : component of normal vector of test solution in y direction
    wavetype: Choose the wave type "P" or "S"
    viz     : Vizualize results if true
    savefile: Save plotfiles if true

    OUTPUT:
    Returns the absolute value of the error in all node points
    """
    # Files for saving
    file1 = File("%s_str-%s-xel-%s-yel-%s-dt-%s-nx-%s-ny-%s.pvd"\
                 % (wavetype, savefile, xel, yel, dt, nx, ny))
    
    # Set solver parameters
    solver = LUSolver("mumps")
    
    # Compute number of elements in x and y direction
    l = L*xel
    m = h*yel
    t = 0

    # Function space and functions
    mesh = RectangleMesh(0,0,L,h,l,m)
    V = VectorFunctionSpace(mesh, "CG", 1)
    Vf = FunctionSpace(mesh, "CG", 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    # Constants
    stepr = Constant(dt**2/rho)

    # Wave type
    if wavetype == "P":
        Au = A*nx
        Av = A*ny
        vel = (lamda + 2*mu)/rho*(nx**2+ny**2)
        k = omega/mt.sqrt(vel)
        g = Expression(("-2*mu*A*k*nx*ny*sin(k*nx*x[0]+k*ny*x[1]-omega*t)",
                        """-lamda*A*k*nx*nx*sin(k*nx*x[0]+k*ny*x[1]-omega*t)
                        -lamda*A*k*ny*ny*sin(k*nx*x[0]+k*ny*x[1]-omega*t)
                        -2*mu*A*k*ny*ny*sin(k*nx*x[0]+k*ny*x[1]-omega*t)"""),
                       mu=mu,A=A,k=k,nx=nx,ny=ny,omega=omega,lamda=lamda,
                       t=t)
                            
        
    elif wavetype == "S":
        Au = A*ny
        Av = -A*nx
        vel = mu/rho*(nx**2+ny**2)
        k = omega/mt.sqrt(vel)
        g = Expression(("""mu*A*k*nx*nx*sin(k*(nx*x[0]+ny*x[1])-omega*t)
                        -mu*A*k*ny*ny*sin(k*(nx*x[0]+ny*x[1])-omega*t)""",
                        "2*mu*A*k*nx*ny*sin(k*(nx*x[0]+ny*x[1])-omega*t)"),
                       mu=mu,A=A,k=k,nx=nx,ny=ny,omega=omega,lamda=lamda,
                       t=t)

    # Initial conditions
    Ixy = Expression(("Au*cos(k*nx*x[0] + k*ny*x[1] - omega*t)",
                      "Av*cos(k*nx*x[0] + k*ny*x[1] - omega*t)"),
                     Au=Au,Av=Av,nx=nx,ny=ny,k=k,omega=omega,t=t)
    
    Vxy = Expression(("Au*cos(k*nx*x[0] + k*ny*x[1] - omega*t)",
                      "Av*cos(k*nx*x[0] + k*ny*x[1] - omega*t)"),
                     Au=Au,Av=Av,nx=nx,ny=ny,k=k,omega=omega,t=t+dt)
    
    u2 = interpolate(Ixy, V)
    u1 = interpolate(Vxy, V)
    
    # Set Dirichlet boundary conditions
    def left(x, on_b): return on_b and abs(x[0]) < DOLFIN_EPS
    def bott(x, on_b): return on_b and abs(x[1]) < DOLFIN_EPS
    def righ(x, on_b): return on_b and abs(x[0] - L) < DOLFIN_EPS

    # Set dirichlet values
    lbc = DirichletBC(V, Ixy, left)
    bbc = DirichletBC(V, Ixy, bott)
    rbc = DirichletBC(V, Ixy, righ)

    # List of dirichlet conditions
    bcs = [rbc, bbc, lbc]
            
    # Stress tensor
    def sigma(v):
        return lamda*div(v)*Identity(2) + \
            mu*(grad(v) + grad(v).T)
    
    # Variational forms
    F = inner(u, v)*dx - 2*inner(u1, v)*dx + inner(u2, v)*dx +\
        stepr*inner(sigma(u1), grad(v))*dx - stepr*dot(g, v)*ds
    
    A = assemble(lhs(F))
    u = Function(V)
    t = 2*dt
    ue = Function(V)
    d = Function(V)
    dxx = d.sub(0)
    dyy = d.sub(1)

    # Main loop
    while t <= T + 2*dt + 0.1*dt + DOLFIN_EPS:
        # Update time dependent bc functions
        Ixy.t = t
        g.t = t-dt
        ue.assign(interpolate(Ixy, V))
        
        # Solve
        begin("Solving at time t=%g" %t)
        b = assemble(rhs(F))
        [bc.apply(A, b) for bc in bcs]
        solver.solve(A, u.vector(), b)
        end()
        d.vector()[:] = abs(ue.vector().array() - u.vector().array())

        # Plot solution
        if viz=='solution':
            plot(u, range_max = 1.0, range_min = -1.0,
                 title="Numerical solution")

        if viz == 'xerror':
            plot(dxx, range_min=0.0, range_max=1e-6, mode='color')

        if viz == 'yerror':
            plot(dyy, range_min=0.0, range_max=1e-6, mode='color')

        if savefile == 'solution':
            file1 << u
        elif savefile == 'xerror':
            file1 << dxx
        elif savefile == 'yerror':
            file1 << dyy
        
        u2.assign(u1)
        u1.assign(u)

        t += dt

    # Exact solution
    Ixy.t = t-dt
    uexact = interpolate(Ixy, V)
    
    # Error at time T
    d.vector()[:] = uexact.vector().array() - u.vector().array()
    error = abs(uexact.vector().array() - u.vector().array())
    return error
    
def run_simulation():
    L = 1
    h = 1
    xel = 24
    yel = 24
    dt = 0.0075
    T = 10
    lamda = 1.
    mu = 1.
    rho = 1.
    A = 1.
    omega = 0.5
    nx = 1
    ny = 0
    wavetype = "S"
    viz = 'none'
    savefile = 'yerror'
    error = solver(L,h,xel,yel,dt,T,\
                       lamda,mu,rho,A,omega,nx,ny,wavetype,viz,\
                       savefile)
    norm = mt.sqrt(sum(error**2/len(error)))
    print 20*'--'
    print 'MAXIMUM ERROR'
    print error.max()
    print 20*'--'
    print 'L2 NORM' 
    print norm
    print 20*'--'

def test_convergence():
    L = 1
    h = 1
    T = 5
    lamda = 1.
    mu = 1.
    rho = 1.
    A = 1.
    omega = 0.5
    nx = 0
    ny = 1
    wavetype = "S"

    dtlist = [0.0075, 0.00375, 0.001875]
    xelist = [24, 48, 96]
    yelist = [24, 48, 96]
    
    # Compute errors
    errorlist = []
    l2normlist = []
    for k in range(len(dtlist)):
        dt = dtlist[k]
        xel = xelist[k]
        yel = yelist[k]
        error = solver(L,h,xel,yel,dt,T,\
                           lamda,mu,rho,A,omega,nx,ny,wavetype,viz=False,\
                           savefile=True)

        # Compute l2 norm
        l2 = mt.sqrt(sum(error**2/len(error)))
        l2normlist.append(l2)
        errorlist.append(error.max())
            

    # Check convergence
    cmax = []
    cl2n = []
    for i in range(len(errorlist)-1):
        cmax.append(errorlist[i+1]/errorlist[i])
        cl2n.append(l2normlist[i+1]/l2normlist[i])

    print 40*'--'
    print 'MAXIMUM ERROR'
    print errorlist
    print 40*'--'
    print 'L2 NORM'
    print l2normlist
    print 40*'--'
    print 'CONVERGENCE MAXIMUM ERROR'
    print cmax
    print 40*'--'
    print 'CONVERGENCE L2 NORM'
    print cl2n
    print 40*'--'



def main():
    run_simulation()
    #test_convergence()
    
if __name__ == '__main__':
    main()
