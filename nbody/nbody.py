"""
This program solves the n-body problem for masses in space. The
program takes the n masses and their initial positions and velocities
and computes their movements in time. The program outputs a data file
with the movements of all the masses.
"""
from scitools.std import *
from body import *

class Nbodysolver:
    """
    class for initiating the nbody solver. class takes the file
    with all initial data, and adds to the list of bodies. class
    then calls a solver method for running a solver algorithm
    """
    def __init__(self):
        self.bodies = []

    def readFile(self, filename):
        """
        Reads bodies from a file and and insert to the bodies list. If the
        bodies list already has elements inside, the new elements are
        added at the end.
        """
        infile = open(str(filename + '.txt'), 'r')
        for line in infile:
            if line.startswith("#"):
                name = line[2:len(line)-1]
                mass = float(((infile.next()).split())[1])
                xpos = float(((infile.next()).split())[1])
                ypos = float(((infile.next()).split())[1])
                xvel = float(((infile.next()).split())[1])
                yvel = float(((infile.next()).split())[1])
                self.add(Body(name, mass, array([xpos,ypos]), \
                              array([xvel,yvel])))
        infile.close()
        
    def add(self, body):
        """
        Adds a body to the end of the list. The method checks for
        inconsistencies in the new body, and terminates if they are
        found.
        """
        self.bodies.append(body)

    def remove(self, param):
        """
        Removes an element from the list. The param can either be an
        integer or a name. The program checks if the param is equal to
        a parameter in the bodies list, and removes the element if
        this is true.
        """
        
    def solve(self, dt, T):
        """
        Solves the nbody problem using the euler cromer method for the
        time step dt and the total time T
        """
        # Number of iterations
        n = int(round(T/dt))

        # Newtons constant
        g = 6.67e-11
        bodies = self.bodies
        t = 0
        while t <= T+0.1*dt:
            # Iterate all bodies
            for i in range(len(bodies)):
                bi = bodies[i]
                F = array([0,0])
                j = 0
                
                # Compute forces on body i from bodies j
                while j <= len(bodies) and j != i:
                    bj = bodies[j]
                    rji = bj.position - bi.position
                    r = linalg.norm(rji)
                    F = F + g*bi.mass*bj.mass*rji/r**3
                
                    j += 1
                    
                    # Compute velocities and positions for body i
                    a = F/bi.mass
                    bi.velocity = bi.velocity + dt*a
                    bi.position = bi.position + dt*bi.velocity
                    
            t += dt

def main():
    filename = "mainfile"
    solver = Nbodysolver()
    solver.readFile(filename)
    solver.solve(0.01, 100)
    

if __name__ == '__main__':
    main()
