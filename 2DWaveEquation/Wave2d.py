from scitools.std import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

class WaveSolver2d:
    """
    Super class for initializing the 2d wave solver. Here, you set
    the initial conditions and the rectangular domain. The solver
    assumes dirichlet boundary conditions and handles this by using
    ghost cells outside the boundaries

    PARAMETERS
    ------------------------------------------------------
    I : func(x,y)  Initial condition
    V : func(x,y)  Initial condition
    l : int        Beginning of domain in x-direction
    h : int        Beginning of domain in y-direction
    L : int        End of domain in x-direction
    H : int        End of domain in y-direction
    n : int        Number of gridpoints in x-direction
    m : int        Number of gridpoints in y-direction
    dt: float      Time step for solver
    T : float      The total simulation time
    q : func(x,y)  Variable wave velocity
    b : float      Damping parameter
    f : func(x,y)  Source term
    ------------------------------------------------------
    """

    def setDomain(self, l, h, L, H, n, m):
        """
        Function sets the domain for the solver
        """
        self.l = l
        self.h = h
        self.L = L
        self.H = H
        self.n = n
        self.m = m

    def setParameters(self, q, b, f):
        """
        Function sets the wave parameters and the source term
        """
        self.q = q;
        self.b = b;
        self.f = f

    def setInitialConditions(self, I, V):
        """
        Function sets the Initial Conditions for the solver
        """
        self.I = I; self.V = V
            
    def solve(self, dt, T):
        """
        Solves the wave equation
        """
        l = self.l; h = self.h; L = self.L; H = self.H
        n = self.n; m = self.m
        q = self.q; b = self.b; f = self.f
        I = self.I; V = self.V
        self.dt = dt; self.T = T
        
        # Set the vectors for the current displacements and
        # the previous displacements -dt and -2dt 
        self.un = zeros((n+3, m+3))
        self.u1 = zeros((n+3, m+3))
        self.u2 = zeros((n+3, m+3))
        self.Ia = zeros((n+3, m+3))
        self.fa = zeros((n+3, m+3))
        self.qa = zeros((n+3, m+3))
        self.Va = zeros((n+3, m+3))

        # Define the grid spacings
        self.dx = float(L)/n
        self.dy = float(L)/m

        # Index sets
        self.Ix = range(1, self.un.shape[0]-1)
        self.Iy = range(1, self.un.shape[1]-1)

        # Define x and y directions
        self.x = linspace(0,L,n+3)
        self.y = linspace(0,H,m+3)
        self.xv = self.x[1:-1, newaxis]
        self.yv = self.y[newaxis, 1:-1]

        
        # Special help values
        self.A = 1./(b*dt/2. + 1)
        self.B = b*dt/2. - 1
        self.dxt = 1./(2.*self.dx**2)
        self.dyt = 1./(2.*self.dy**2)

        # Set function values
        self.qa[1:-1,1:-1] = q(self.xv,self.yv)
        self.Ia[1:-1,1:-1] = I(self.xv,self.yv)
        self.Va[1:-1,1:-1] = V(self.xv,self.yv)
        self.fa[1:-1,1:-1] = f(self.xv,self.yv,0)
        
        # Ghost cells for q
        self.qa[1:-1,0] = self.qa[1:-1,2]
        self.qa[1:-1,m+2] = self.qa[1:-1,m]
        self.qa[0,1:-1] = self.qa[2,1:-1]
        self.qa[n+2,1:-1] = self.qa[n,1:-1]

        # Check stability condition
        self.checkStability()
        
        # Apply initial conditions
        self.t = 2*dt
        u2, u1 = self.applyInitialConditions()

        # Construct figure and grid
        plt.ion()
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        X,Y = meshgrid(self.x[1:-1], self.y[1:-1])
        
        # Start loop
        while self.t < T + 0.1*dt:

            # Update u
            self.un = self.advance()
            self.t += dt

        raw_input("Enter something to exit: ")

    def checkStability(self):
        """
        Checks if the stability condition for the wave equation is satisfied
        """
        if self.dt <= 1./sqrt(abs(self.qa.max()))*\
           1./sqrt(1./self.dx**2 + 1./self.dy**2):
            print "####################################################"
            print "WARNING!!!"
            print "Stability condition is not satisfied."
            print "Implementation might fail"
            print "####################################################"
        else:
            print "####################################################"
            print "Stability condition is satisfied"
            print "####################################################"
            
        
class VectorizedWaveSolver(WaveSolver2d):
    """
    Class solves the wave equation using vectorized computations.
    The method applyInitialConditions sets the initial conditions to
    the domain. The method solve
    """
    def applyInitialConditions(self):
        """
        Apply initial conditions to domain using vectorized computations
        """
        # First initial condition
        self.u2[1:-1,1:-1] = self.Ia[1:-1,1:-1]

        # Update ghost points
        self.u2[1:-1,0] = self.u2[1:-1,2]
        self.u2[1:-1,self.m+2] = self.u2[1:-1,self.m]
        self.u2[0,1:-1] = self.u2[2,1:-1]
        self.u2[self.n+2,1:-1] = self.u2[self.n,1:-1]
        
        # Second initial condition
        self.fa[1:-1,1:-1] = self.f(self.x[1:-1], self.y[1:-1], self.t+self.dt)
            
        qux = self.dxt*((self.qa[1:-1,1:-1] + self.qa[2:,1:-1])*\
                   (self.u2[2:,1:-1] - self.u2[1:-1,1:-1]) - \
                   (self.qa[:-2,1:-1] + self.qa[:-2,1:-1])*\
                   (self.u2[1:-1,1:-1] - self.u2[:-2,1:-1]))
        
        quy = self.dyt*((self.qa[1:-1,1:-1] + self.qa[1:-1,2:])*\
                   (self.u2[1:-1,2:] - self.u2[1:-1,:-2]) - \
                   (self.qa[1:-1,1:-1] + self.qa[1:-1,:-2])*\
                   (self.u2[1:-1,1:-1] - self.u2[1:-1,:-2]))
        
        self.u1[1:-1,1:-1] = self.u2[1:-1,1:-1] + self.dt*self.Va[1:-1,1:-1] + \
                        self.b*self.dt**2/2.*self.Va[1:-1,1:-1] + \
                        1/2.*self.dt**2*(qux + quy + self.fa[1:-1,1:-1])
            
        # Update ghost values
        self.u1[1:-1,0] = self.u1[1:-1,2]
        self.u1[1:-1,self.m+2] = self.u1[1:-1,self.m]
        self.u1[0,1:-1] = self.u1[2,1:-1]
        self.u1[self.n+2,1:-1] = self.u1[self.n,1:-1]

        return self.u2, self.u1
        
    def advance(self):
        """
        Advance the solution to the next step
        """
        self.fa[1:-1,1:-1] = self.f(self.x[1:-1], self.y[1:-1], self.t)
                               
        qux = self.dxt*((self.qa[1:-1,1:-1] + self.qa[2:,1:-1])*\
                   (self.u1[2:,1:-1] - self.u1[1:-1,1:-1]) - \
                   (self.qa[:-2,1:-1] + self.qa[:-2,1:-1])* \
                   (self.u1[1:-1,1:-1] - self.u1[:-2,1:-1]))
                           
        quy = self.dyt*((self.qa[1:-1,1:-1] + self.qa[1:-1,2:])* \
                   (self.u1[1:-1,2:] - self.u1[1:-1,1:-1]) - \
                   (self.qa[1:-1,1:-1] + self.qa[1:-1,:-2])* \
                   (self.u1[1:-1,1:-1] - self.u1[1:-1,:-2]))
        
        self.un[1:-1,1:-1] = self.A*(2*self.u1[1:-1,1:-1] + \
                                     self.u2[1:-1,1:-1]*self.B + \
                                     self.dt**2*(qux + quy + \
                                                 self.fa[1:-1,1:-1]))

        # Update ghost values
        self.un[1:-1,0] = self.un[1:-1,2]
        self.un[1:-1,self.m+2] = self.un[1:-1,self.m]
        self.un[0,1:-1] = self.un[2,1:-1]
        self.un[self.n+2,1:-1] = self.un[self.n,1:-1]

        # Update values for next iteration
        self.u2[:,:] = self.u1[:,:]
        self.u1[:,:] = self.un[:,:]

        return self.un

def main():
    solver = VectorizedWaveSolver()
    def I(x,y): return exp(-0.5*((x - 0.5))**2 - 0.5*((y - 0.5))**2)
    def V(x,y): return 0
    def q(x,y): return 1
    def f(x,y,t): return 0
    
    solver.setDomain(0, 0, 1, 1, 101, 101)
    solver.setParameters(q, 1, f)
    solver.setInitialConditions(I, V)
    solver.solve(0.01, 10)

if __name__ == '__main__':
    main()
