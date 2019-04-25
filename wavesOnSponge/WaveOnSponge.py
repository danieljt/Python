from dolfin import *
import math as mt

def solver(L,h,xel,yel,xs,dt,T,omega,vel,k,damp,viz,save):
    """
    Function for solving the scalar wave equation in a
    rectangular domain with given boundary and initial
    conditions
    -------------------------------------------------
    INPUT:
    L:     Length of domain
    h:     Height of domain
    xel:   Number of elements per unit length in the x-direction
    yel:   Number of elements per unit length in the y-direction
    xs:    Coordinate of the vertical line seperating fluid and sponge
    dt:    Time step
    T:     Total simulation time
    omega: Angular frequency
    vel:   Velocity of waves
    k:     Constant determinging the
    damp:  lin for linear, and quad for quadratic damping
    viz:   True for showing simulation plot
    save:  True for saving errors at time T
    --------------------------------------------------
    OUTPUT:
    Returns the error between analytic and exact
    solution in the fluid layer
    Saves plots of the component errors if save=True
    """
    # Create file for storing file
    savefile = File("viz-%s-%s-%s-%s-%s-%s-%s.pvd" \
                    % (damp,L,h,xel,xs,dt,T))
    
    # Starting time
    t = 0
    
    # elements per length
    l = L*xel
    m = h*yel

    # Define functionspace
    mesh = RectangleMesh(0,0,L,h,l,m)
    V = FunctionSpace(mesh, "CG", 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    # Define subdomains
    class Fluid(SubDomain):
        def inside(self,x,on_boundary):
            return (between(x[0], (0,xs)))

    class Sponge(SubDomain):
        def inside(self,x,on_boundary):
            return (between(x[0], (xs,L)))

    fluid = Fluid()
    sponge = Sponge()                    
    domains = CellFunction("size_t", mesh)
    domains.set_all(0)
    fluid.mark(domains, 0)
    sponge.mark(domains, 1)

    # Create submesh from fluid domain
    submesh = SubMesh(mesh, fluid)
    Vf = FunctionSpace(submesh, "CG", 1)
    ud = Function(Vf)
    ue = Function(Vf)
    uf = Function(Vf)
    
    # Variable expressions
    ce = Constant(vel)

    # Set the damping to lin or quad
    if damp=="lin":
        be = Expression("x[0] < xs ? 0 : 10*(x[0]-xs)", xs=xs)
    elif damp=="quad":
        be = Expression("x[0] < xs ? 0 : 10*(pow(x[0],2)-2*xs*x[0]+pow(xs,2))",
                        xs=xs)
    else:
        print "Insert lin or quad"
        exit()
    
    # Define important constants
    step2 = Constant(1/dt**2)
    step3 = Constant(1/(2*dt))

    # Initial conditions
    Ixy = Constant(0)
    Vxy = Constant(0)
    
    # Essential boundary conditions
    inflow = Expression("sin(omega*t)*cos(pi*x[1]/(2*h)*(1 + k))",
                        omega=omega,h=h,k=k,t=t)
    free = Constant(0)
    def surface(x, ob): return ob and abs(x[1]-h) < DOLFIN_EPS
    def leftfun(x, ob): return ob and abs(x[0]) < DOLFIN_EPS

    left = DirichletBC(V, inflow, leftfun)
    topp = DirichletBC(V, free, surface)
    bcs = [left, topp]

    # Set all functions into domain
    c = interpolate(ce, V)
    b = interpolate(be, V)
    u1 = interpolate(Ixy, V)
    u2 = interpolate(Vxy, V)

    # Exact solution
    lk = mt.pi*(1 + k)/(2.*h)
    kk = mt.sqrt(omega**2/vel**2 - lk**2)
    uexact = Expression("sin(omega*t - kk*x[0])*cos(lk*x[1])",
                        omega=omega,kk=kk,lk=lk,t=t-2*dt)
    
    # Variational forms
    F = step2*inner(u,v)*dx - 2*step2*inner(u1,v)*dx + step2*inner(u2,v)*dx +\
        b*step3*inner(u,v)*dx - b*step3*inner(u2,v)*dx +\
        c*inner(nabla_grad(u1), nabla_grad(v))*dx

    A = assemble(lhs(F))
    u = Function(V)
    t = 2*dt
    
    while t <= T + 2*dt + DOLFIN_EPS:
        inflow.t = t
        uexact.t = t-2*dt
        begin("Computing at time level t = %g" %t)
        LL = assemble(rhs(F))
        [bc.apply(A,LL) for bc in bcs]
        solve(A, u.vector(), LL)
        end()
        ue.assign(interpolate(uexact, Vf))
        uf.assign(interpolate(u2, Vf))
        ud.vector()[:] = abs(ue.vector().array() - uf.vector().array())
                  
        if viz == 'solution':
            plot(u2, range_max=1.0, range_min=-1.0, title="Numerical solution")
            if save == True:
                savefile << u2
                  
        if viz == 'error':
            plot(ud, rescale=False, mode='color', title='Error in fluid layer')
            if save == True:
                savefile << ud
            
        u2.assign(u1)
        u1.assign(u)

        t += dt

    return ud.vector().array()
    

def run_simulation():
    """
    Test program for running an experiment showing the
    plot on screen with given values and returning the
    error. The maximum and L2 norm errors
    are printed at terminal
    """
    L = 3
    h = 1
    xel = 24
    yel = 24
    xs = 1
    dt = 0.01
    T = 10
    omega = 10.
    vel = 1.
    k = 0
    damp="lin"
    viz = 'error'
    save = True
    error = solver(L,h,xel,yel,xs,dt,T,omega,vel,k,damp,viz,save)
    error_max = error.max()
    error_l2n = mt.sqrt(sum(error**2/len(error)))

    print "Maximum error: ", error_max
    print "L2 norm error: ", error_l2n

def test_convergence():
    """
    Program for running a convergence test with given physical
    values. The time and spatial steps are halved to test that
    convergence is reached. Component errors are then saved to VTK
    files.
    """
    L = 3
    h = 1
    xs = 1
    T = 10
    vel = 1.
    omega = 10.
    k = 0
    damp = "quad"
    viz=False
    save=True

    # Lists to store error values
    E_max = []
    E_l2n = []

    # Lists with dt, dx and dy values
    timestep = [0.01, 0.005, 0.0025]
    xelement = [24, 48, 96]
    yelement = [24, 48, 96]
    
    for i in range(len(timestep)):
        dt = timestep[i]
        xel = xelement[i]
        yel = yelement[i]
        error = solver(L,h,xel,yel,xs,dt,T,omega,vel,k,damp,viz,save)

        error_max = error.max()
        error_l2n = mt.sqrt(sum(error**2/len(error)))

        E_max.append(error_max)
        E_l2n.append(error_l2n)
        

    # Check convergence
    C_max = []
    C_l2n = []
    for i in range(len(E_max)-1):
        C_max.append(E_max[i+1]/E_max[i])
        C_l2n.append(E_l2n[i+1]/E_l2n[i])
        
    print 40*'--'
    print 'MAXIMUM ERROR'
    print E_max
    print 40*'--'
    print 'L2 NORM'
    print E_l2n
    print 40*'--'
    print 'CONVERGENCE MAXIMUM ERROR'
    print C_max
    print 40*'--'
    print 'CONVERGENCE L2 NORM'
    print C_l2n
    print 40*'--'

        

    
    
def main():
    #run_simulation()
    test_convergence()

if __name__=='__main__':
    main()
