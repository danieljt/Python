"""
This program aims at computing the movement of diffent bodies in 2D space by
using newtons laws of motion. For simplicity, relativistic effects are
neglected and the bodies are assumed circular. The program is divided
into 3 classes.

Body   : Holds information about the bodies
Solver : Holds the solver
Nbody  : Main class for 
"""

from scitools.std import *

class Body:
    def __init__(self, mass, radius, position, velocity):
        self.m = mass
        self.r = radius
        self.p = position
        self.v = velocity
