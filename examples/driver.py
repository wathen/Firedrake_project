from firedrake import *
from project import Poisson, ExactSol, ConvecDiff, Helmholtz, misc
import numpy as np
from matplotlib.pylab import plt


for i in range(4):
    n = 2**(i + 2)
    mesh = UnitCubeMesh(n, n, n)
    V = VectorFunctionSpace(mesh, 'CG', 1)
    p = Helmholtz(V, nu=1., opt={'form': 'linear'})
    p.uE = ['sin(y)*cos(x)', 'y*x', 'tan(z)*y*x']
    p.k = 1.0
    solve(p.form - p.rhs == 0, p.sol, bcs=p.bcs)
    print(errornorm(misc.Expression(p.uE, V, 1), p.sol))
    e = Function(V)
    e.dat.data[:] = p.sol.dat.data - misc.Expression(p.uE, V).dat.data
