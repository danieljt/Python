from scitools.std import *
import matplotlib.pyplot as plt

def solver(nx,ny,L,H,dt,T,Re,M,U):
    """
    Function for evaluating the 2D lid driven cavity flow from the
    maccormack scheme in a rectangular domain. The solver is
    non-dimensionalized and assumes a constant grid spacing in
    the x- and y directions. The solver uses vectorized code

    INPUT:
    nx: Number of gridpoints in x-direction
    ny: Number of gridpoints in y-direction
    L : The Length of the domain
    H : The Height of the domain
    dt: The time step
    T : The total simulation time
    Re: The Reynolds number
    M : Mach number
    U : Velocity of lid
    """
    # Domain arrays
    dx = float(L)/nx
    dy = float(H)/ny
    x, y = mgrid[0:L+dx:dx, 0:H+dy:dy]
    
    # time
    N = int(round(T/float(dt)))
    t = linspace(0,N*dt,N+1)
    
    # Constants 
    a1 = dt/dx
    a2 = dt/dy
    a3 = dt/(dx*M**2)
    a4 = dt/(dy*M**2)
    a5 = 4*dt/(3*Re*dx**2)
    a6 = dt/(Re*dy**2)
    a7 = dt/(Re*dx**2)
    a8 = 4*dt/(3*Re*dy**2)
    a9 = dt/(12*Re*dx*dy)
    a10 = 2*(a5 + a6)
    a11 = 2*(a7 + a8)

    # Matrixes for values
    r = zeros((nx+1,ny+1))
    rs = zeros((nx+1,ny+1))
    rn = zeros((nx+1,ny+1))
    u = zeros((nx+1,ny+1))
    us = zeros((nx+1,ny+1))
    un = zeros((nx+1,ny+1))
    v = zeros((nx+1,ny+1))
    vs = zeros((nx+1,ny+1))
    vn = zeros((nx+1,ny+1))
    ur = zeros((nx+1,ny+1))
    vr = zeros((nx+1,ny+1))
    urs = zeros((nx+1,ny+1))
    vrs = zeros((nx+1,ny+1))
    urn = zeros((nx+1,ny+1))
    vrn = zeros((nx+1,ny+1))

    # Set initial conditions
    r[:,:] = 1
    ur[:,-1] = U*r[:,-1]

    # Start the solver loop
    for i in range(len(t)):
        # Calculate velocities at current time step
        u[:,:] = ur[:,:]/r[:,:]
        v[:,:] = vr[:,:]/r[:,:]

        # Predictor step
        rs,urs,vrs = predictor(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,\
                               r,rs,u,ur,urs,v,vr,vrs,U)
        
        # Calculate predicted velocities
        us = urs[:,:]/rs[:,:]
        vs = vrs[:,:]/rs[:,:]

        # Corrector step
        rn,urn,vrn =  corrector(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,\
                                r,rs,rn,ur,us,un,urs,urn,
                                vr,vs,vn,vrs,vrn,U)

        # update values for next time step
        r[:,:] = rn[:,:]
        ur[:,:] = urn[:,:]
        vr[:,:] = vrn[:,:]

    # Create vector plot
    figure = plt.figure()
    plt.xlabel("x")
    plt.ylabel("y")
    plt.quiver(x,y,u,v)
    plt.axis([0-5*dx,L+5*dx,0-5*dy,H+5*dy])
    raw_input("Enter something to exit: ")
    
def predictor(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,\
              r,rs,u,ur,urs,v,vr,vrs,U):
    """
    This function computes the predictor step 
    """
    # Interior nodes
    rs[1:-1,1:-1] = r[1:-1,1:-1] - \
                    a1*(ur[2:,1:-1] - ur[1:-1,1:-1]) -\
                    a2*(vr[1:-1,2:] - vr[1:-1,1:-1])
    
    urs[1:-1,1:-1] = ur[1:-1,1:-1] - \
                     a3*(r[2:,1:-1] - r[1:-1,1:-1]) -\
                     a1*(ur[2:,1:-1]*u[2:,1:-1] - \
                         ur[1:-1,1:-1]*u[1:-1,1:-1]) -\
                     a2*(ur[1:-1,2:]*v[1:-1,2:] - \
                         ur[1:-1,1:-1]*v[1:-1,1:-1]) -\
                     a10*u[1:-1,1:-1] + a5*(u[2:,1:-1] + u[:-2,1:-1]) +\
                     a6*(u[1:-1,2:] + u[1:-1,:-2]) +\
                     a9*(v[2:,2:] + v[:-2,:-2] - v[2:,:-2] - v[:-2,2:])

    
    vrs[1:-1,1:-1] = vr[1:-1,1:-1] - a4*(r[1:-1,2:] - r[1:-1,1:-1]) -\
                     a1*(ur[2:,1:-1]*v[2:,1:-1] - \
                         ur[1:-1,1:-1]*v[1:-1,1:-1]) -\
                     a2*(vr[1:-1,2:]*v[1:-1,2:] -\
                         vr[1:-1,1:-1]*v[1:-1,1:-1]) -\
                     a11*v[1:-1,1:-1] + a7*(v[2:,1:-1] + v[:-2,1:-1]) +\
                     a8*(v[1:-1,2:] + v[1:-1,:-2]) +\
                     a9*(u[2:,2:] + u[:-2,:-2] - u[2:,:-2] - u[:-2,2:])

    # Update Boundary values for density
    rs[0,:] = r[0,:] - 0.5*a1*(-ur[2,:]+4*ur[1,:]-3*ur[0,:])
    rs[-1,:] = r[-1,:] - 0.5*a1*(-ur[-3,:]+4*ur[-2,:]-3*ur[-1,:])
    rs[:,0] = r[:,0] - 0.5*a2*(-vr[:,2] + 4*vr[:,1] -3*vr[:,0])
    rs[1:-1,-1] = r[1:-1,-1] - 0.5*U*a1*(r[2:,-1] - r[:-2,-1]) +\
               0.5*a2*(-vr[1:-1,-3] + 4*vr[1:-1,-2] -3*vr[1:-1,-1])

    # Update Boundary values for velocity components
    vrs[0,:] = 0; vrs[-1,:] = 0; vrs[:,0] = 0; vrs[:,-1] = 0
    urs[0,:] = 0; urs[-1,:] = 0; urs[:,0] = 0; urs[:,-1] = U*rs[:,-1]
    
    return rs, urs, vrs
    

def corrector(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,\
              r,rs,rn,ur,us,un,urs,urn,
              vr,vs,vn,vrs,vrn,U):
    """
    This function computes the corrector step
    """
    # Update Interior nodes
    rn[1:-1,1:-1] = 0.5*(r[1:-1,1:-1] + rs[1:-1,1:-1]) -\
                    0.5*a1*(urs[1:-1,1:-1] - urs[:-2,1:-1]) -\
                    0.5*a2*(vrs[1:-1,1:-1] - vrs[1:-1,:-2])

    urn[1:-1,1:-1] = 0.5*(ur[1:-1,1:-1] + urs[1:-1,1:-1]) -\
                     0.5*a3*(rs[1:-1,1:-1] - rs[:-2,1:-1]) -\
                     0.5*a1*(urs[1:-1,1:-1]*us[1:-1,1:-1] -\
                             urs[:-2,1:-1]*us[:-2,1:-1]) -\
                     0.5*a2*(urs[1:-1,1:-1]*us[1:-1,1:-1] -\
                             urs[1:-1,:-2]*us[1:-1,:-2]) -\
                     0.5*a10*us[1:-1,1:-1] +\
                     0.5*a5*(us[2:,1:-1] + us[:-2,1:-1]) +\
                     0.5*a6*(us[1:-1,2:] + us[1:-1,:-2]) +\
                     0.5*a9*(vs[2:,2:] + vs[:-2,:-2] -\
                             vs[2:,:-2] - vs[:-2,2:])

    vrn[1:-1,1:-1] = 0.5*(vr[1:-1,1:-1] + vrs[1:-1,1:-1]) -\
                     0.5*a4*(rs[1:-1,1:-1] - rs[1:-1,:-2]) -\
                     0.5*a1*(urs[1:-1,1:-1]*vs[1:-1,1:-1] -\
                             urs[:-2,1:-1]*vs[:-2,1:-1]) -\
                     0.5*a2*(vrs[1:-1,1:-1]*vs[1:-1,1:-1] -\
                             vrs[1:-1,:-2]*vs[1:-1,:-2]) -\
                     0.5*a11*vs[1:-1,1:-1] +\
                     0.5*a7*(vs[2:,1:-1] + vs[:-2,1:-1]) +\
                     0.5*a8*(vs[1:-1,2:] + vs[1:-1,:-2]) +\
                     0.5*a9*(us[2:,2:] + us[:-2,:-2] -\
                             us[2:,:-2] - us[:-2,2:])

    # Update Boundary conditions for density
    rn[0,:] = 0.5*(r[0,:] + rs[0,:]) -\
              0.25*a1*(-urs[2,:] + 4*urs[1,:] - 3*urs[0,:])
    rn[-1,:] = 0.5*(r[0,:] + rs[0,:]) + \
               0.25*a1*(-urs[-3,:] + 4*urs[-2,:] - 3*urs[-1,:])
    rn[:,0] = 0.5*(r[:,0] + rs[:,0]) -\
              0.25*a2*(-vrs[:,2] + 4*vrs[:,1] - 3*vrs[:,0])
    rn[1:-1,-1] = 0.5*(r[2:,-1] + rs[:-2,-1]) -\
               0.25*U*a1*(rs[2:,-1] - rs[:-2,-1]) +\
               0.25*a2*(-vrs[1:-1,-3] + 4*vrs[1:-1,-2] - 3*vrs[1:-1,-1])

    # Update boundary conditions for velocity components
    vrn[0,:] = 0; vrn[-1,:] = 0; vrn[:,0] = 0; vrn[:,-1] = 0
    urn[0,:] = 0; urn[-1,:] = 0; urn[:,0] = 0; urn[:,-1] = U*rn[:,-1]

    return rn, urn, vrn

def test():
    # Test solver
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

if __name__ == '__main__':
    test()
