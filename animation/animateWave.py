"""
Test program for animating wavemotion using
matplotlibs animation module. 
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def f(x,t):
    return np.cos(x - 6*t)

class UpdateFunction(object):
    def __init__(self, ax):
        self.line, = ax.plot([],[],'k-')
        self.x = np.linspace(0,5,1001)
        self.ax = ax

        self.ax.set_xlim(0,5)
        self.ax.set_ylim(-2,2)
        self.ax.grid(True)

    def init(self):
        self.line.set_data([],[])
        return  self.line,

    def __call__(self, i):
        if i==0:
            return self.init()

        y = f(self.x,i)
        self.line.set_data(self.x,y)
        return self.line,

fig, ax = plt.subplots()
uf = UpdateFunction(ax)
anim = FuncAnimation(fig, uf, frames = np.arange(100), init_func=uf.init,
                     interval = 100, blit=True)
plt.show()
