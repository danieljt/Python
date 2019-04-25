from dolfin import *
import math as mt

def solver(L,h,xel,yel,dt,T,lamda,mu,rho,A,omega,nx,ny,wavetype,viz,savefile):
    """
    Function for solving the seismic wave equation on a rectangular
    domain with with inhomogeneous dirichlet bcs on all sides. The solver
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
    file1 = File("%s_%s_xyel_%s_yel_%s_dt_%s_nx_%s_ny_%s.pvd" \
                 % (wavetype, savefile, xel, yel, dt, nx, ny))
  
    # Set solver and plotter
    solver = LUSolver("mumps")
    
    # Compute number of elements in x and y direction
    l = L*xel
    m = h*yel

    # Function space and functions
    mesh = RectangleMesh(0,0,L,h,l,m)
    V = VectorFunctionSpace(mesh, "CG", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    ue = Function(V)
    dxy = Function(V)
    dxx = dxy.sub(0)
    dyy = dxy.sub(1)
    
    # Constants
    stepr = Constant(dt**2/rho)

    # Test Wave type
    if wavetype == "P": # Pressure wave
        Au = A*nx
        Av = A*ny
        vel = mt.sqrt((lamda + 2*mu)/rho*(nx**2+ny**2)) # Wave velocity
        k = omega/vel # Dispersion relation
        
    elif wavetype == "S": # Shear wave
        Au = A*ny
        Av = -A*nx
        vel = mt.sqrt(mu/rho*(nx**2 + ny**2)) # Wave velocity
        k = omega/vel # Dispersion relation

    t = 0
    # Initial conditions
    Ixy = Expression(("Au*cos(k*nx*x[0] + k*ny*x[1] - omega*t)",
                      "Av*cos(k*nx*x[0] + k*ny*x[1] - omega*t)"),
                     Au=Au,Av=Av,nx=nx,ny=ny,k=k,omega=omega,t=t)
   
    Vxy = Expression(("Au*cos(k*nx*x[0] + k*ny*x[1] - omega*t)",
                      "Av*cos(k*nx*x[0] + k*ny*x[1] - omega*t)"),
                     Au=Au,Av=Av,nx=nx,ny=ny,k=k,omega=omega,t=t+dt)

    u2 = interpolate(Ixy, V)
    u1 = interpolate(Vxy, V)
        
    # Boundary condition
    def boundary(x, on_boundary): return on_boundary
    bc = DirichletBC(V, Ixy, boundary)
    
    # Stress tensor
    def sigma(u, lamda, mu):
        return lamda*div(u)*Identity(2) + mu*(grad(u) + grad(u).T)

    # Variational form
    F = inner(u, v)*dx - 2*inner(u1, v)*dx + inner(u2, v)*dx +\
        stepr*inner(sigma(u1, lamda, mu), grad(v))*dx

    A = assemble(lhs(F))  # Assemble left hand side
    u = Function(V)
    t = 2*dt

    while t <= T + DOLFIN_EPS:
        # Update time dependent bc functions
        Ixy.t = t
        ue.assign(interpolate(Ixy, V))
        
        # Solve
        begin("Solving at time step t=%g" % t)
        b = assemble(rhs(F))
        bc.apply(A, b)
        solver.solve(A, u.vector(), b)
        end()
        dxy.vector()[:] = abs(ue.vector().array() - u.vector().array())
        
        # Plot solution
        if viz=='solution':
            plot(u, range_max = 1.0, range_min = -1.0,
                 title="Numerical solution")

        elif viz == 'xerror':
            plot(dxx, range_max=1e-6, range_min=-1e-6, mode='color')

        elif viz == 'yerror':
            plot(dyy, range_max=1e-6, range_min=-1e-6, mode='color')

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

    # Compute component differences
    dxy.vector()[:] = uexact.vector().array() - u.vector().array()

    # return the error
    error = abs(uexact.vector().array() - u.vector().array())
    return error    

def run_simulation():
    L = 1
    h = 1
    xel = 24
    yel = 24
    dt = 0.0075
    T = 10.0
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
    error = solver(L,h,xel,yel,dt,T,
                   lamda,mu,rho,A,omega,nx,ny,wavetype,viz,
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
    normlist = []
    for k in range(len(dtlist)):
        dt = dtlist[k]
        xel = xelist[k]
        yel = yelist[k]
        error = solver(L,h,xel,yel,dt,T,
                       lamda,mu,rho,A,omega,nx,ny,wavetype,viz=False,
                       savefile = True)

        # Compute l2 norm
        norm = mt.sqrt(sum((error)**2/(len(error))))
        normlist.append(norm)
        errorlist.append(error.max())
            

    # Check convergence
    cmax = []
    cl2n = []
    for i in range(len(errorlist)-1):
        cmax.append(errorlist[i+1]/errorlist[i])
        cl2n.append(normlist[i+1]/normlist[i])

    print 40*'--'
    print 'MAXIMUM ERROR'
    print errorlist
    print 40*'--'
    print 'L2 NORM'
    print normlist
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
