The nbody programs aim at solving a system of bodies in space driven by gravity.
The program can be used for star systems, solar systems and other
celestial movents in space, and further extensions of the program should be
easy.

-----------------------------------------------------------------
Classes:
-----------------------------------------------------------------
Body    : Contains information about a body
- Mass
- Position
- Velocity

Nbody   :
- Bodylist
- readfile
- solve
- plot

-----------------------------------------------------------------
Usage:
-----------------------------------------------------------------
From Nbody import *

Nbody.readfile("Filename")                     # Add bodies from file
Nbody.addBody(name, mass, position, velocity)  # Add body manually
Nbody.solve(dt, T, solverType)                 # Solve problem
Nbody.plot()                                   # Plot end solution


