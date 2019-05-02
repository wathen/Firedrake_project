from firedrake import *
from project import Poisson, ConvecDiff, Helmholtz, Stokes, ExactSol, misc, FiniteElementSpaces
import numpy as np
from matplotlib.pylab import plt


for i in range(4):
    n = 2**(i + 2)
    mesh = UnitSquareMesh(n, n)
    V = FiniteElementSpaces.TaylorHood(mesh)
    # V = VectorFunctionSpace(mesh, 'CG', 1)
    p = Stokes(V, nu=1., opt={'form': 'bilinear'})
