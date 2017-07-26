The program wave-sponge tries to solve the wave equation with a damping layer.
When modelling propagating values in a very large or infinite material using the finite element method, we have two fundemental problems.

1) It is impossible to implement an infinite material, and it
   it is very expensive to implement a large one. Even if a large layer
   is implemented, the values reaching the end will induce reflections
   that will propagate back and in time destroy the solution

2) The finite element method needs set boundary conditions to be able to 
   solve it's problems. This is very difficult for time dependent problems
   and most likely impossible as it is not possible to predict the behaviour

A solution to this problem is the damping layer, or sponge layer as it's often called. The sponge layer is a small domain outside our main domain which has the purpose of destroying the propagating value. This eliminates the need of implementing very large domains, and prohibits a large portion of the reflected waves that would occur. It's important to note that reflected waves will still be produced. This program illustrates an implementation of such a layer for incoming surface waves with an analytical solution. The program uses P1 elements for the spatial domain and finite differences with second order errors. This should give fourth order convergence rate. So for every time the domains are halved, the error should be reduce to 1/4
