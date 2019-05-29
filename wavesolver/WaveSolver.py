"""
Problem
in the animation part, the initial conditions are not set because
they are only defined in the solve function. Must generalize the
time constants and initial conditions into separate functions.
The animation module seems to work.
"""
from numpy import linspace, zeros, array, asarray, pi, cos, sin, sqrt, log
from matplotlib.pyplot import plot, show, figure
import matplotlib.animation as animation

class Wavesolver:
    def __init__(self, c):
        self.c  = c

    def set_domain(self, x0,x1,nx):
        self.x  = linspace(x0,x1,nx+1)
        self.dx = float((x1 - x0)/(nx))
        self.u0 = zeros(nx+1)
        self.u1 = zeros(nx+1)
        self.u2 = zeros(nx+1)
        self.us = zeros(nx+1)

    def set_time(self, dt, T):
        self.dt = dt
        self.T  = T
        self.n = int(round(T/dt))
        self.t = 0
        
        
    def set_boundary_conditions(self, bcs):
        self.bcs = bcs

    def set_initial_conditions(self, f, g):
        self.f = f
        self.g = g

        self.u0[:] = self.f(self.x[:])
        self.us[:] = self.u0[:]
        self.u1[1:-1] = self.u0[1:-1] + self.dt*self.g(self.x[1:-1]) + \
                        0.5*(self.c*self.dt/self.dx)**2*\
                        (self.u0[2:] - 2*self.u0[1:-1] + self.u0[:-2])
        self.u1[0]  = 0
        self.u1[-1] = 0
        self.t = 2*self.dt
        
    def solve(self):
        while self.t < self.T + 2*self.dt + 0.01*self.dt:
            self.u2= self.advance()

        return self.us

    def advance(self):
        self.u2[1:-1] = 2*self.u1[1:-1] - self.u0[1:-1] + \
                    (self.c*self.dt/self.dx)**2*\
                    (self.u1[2:] - 2*self.u1[1:-1] + self.u1[:-2])

        self.u2[0]  = 0
        self.u2[-1] = 0

        self.u0[:] = self.u1[:]
        self.us[:] = self.u0[:]
        self.u1[:] = self.u2[:]
        self.t += self.dt
        
        return self.u2

    def animate_solution(self, dt, T):
        fig = figure()
        ax = fig.add_subplot(111, autoscale_on=False,
                             xlim=(-self.x[0]-1, self.x[-1]+1),
                             ylim=(-2, 2))
        self.line, = ax.plot([], [])
        movie = animation.FuncAnimation(fig, self.animate,
                                        frames=self.n, blit=True)
        show()
    
    def animate(self, i):
        self.advance()
        self.line.set_data(self.x, self.us[:])
        return self.line,
        
def test_convergence():
    x0 = 0
    x1 = 1
    c  = 1
    T  = 10
    def f(x): return sin(pi*x)
    def g(x): return 0
    def e(x,t,c): return cos(c*pi*t)*sin(pi*x)
    nl = [40, 80, 160, 320]
    dt = [0.025, 0.0125, 0.00625, 0.003125]
    el = []

    for i in range(len(nl)):
        x = linspace(x0,x1,nl[i]+1)
        solver = Wavesolver(c)
        solver.set_domain(x0,x1,nl[i])
        solver.set_time(dt[i], T )
        solver.set_initial_conditions(f, g)
        un = solver.solve()
        el.append((abs(e(x,T,c) - un)).max())

    nl = asarray([20, 40, 80, 160])
    nl = 1./nl
    el = asarray(el)
    r = log(el[:-1]/el[1:])/(log(nl[:-1]/nl[1:]))
    print(el)
    print(r)

def test_animation():
    x0 = 0
    x1 = 1
    c  = 1
    T  = 1
    nl = 20
    dt = 0.05
    def f(x): return sin(pi*x)
    def g(x): return 0
    def e(x,t,c): return cos(c*pi*t)*sin(pi*x)
    solver = Wavesolver(c)
    solver.set_domain(x0,x1,nl)
    solver.set_time(dt, T)
    solver.set_initial_conditions(f, g)
    solver.animate_solution(dt, T)
    #x = linspace(x0,x1,nl+1)
    #plot(x, e(x,T,c))
    show()
    
if __name__ == '__main__':
    #test_convergence()
    test_animation()

    
    
    
