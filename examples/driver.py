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
    p.uE = ['x', 'y']
    p.pE = ['x']
    # p.k = 1.0
    assemble(p.form)
    solve(p.form == p.rhs, p.sol, bcs=p.bcs,
          solver_parameters={'ksp_type': 'preonly',
                             'pc_type': 'lu'})
    print(errornorm(misc.Expression(p.uE, V, 1), p.sol))
    e = Function(V)
    e.dat.data[:] = p.sol.dat.data - misc.Expression(p.uE, V).dat.data
WcmQLSGd2vMWAT11@PH.RC
