from firedrake import *
from project import Poisson
from project import ConvecDiff
from project.utils import Expression
import numpy as np


n = 2**3.
mesh = UnitSquareMesh(n, n)
V = VectorFunctionSpace(mesh, 'CG', 1)

p = ConvecDiff(V, nu=0.01, opt={'form': 'linear'})
p.uE = ['x*x*sin(x)', 'y*x']
p.wind = ['x*x*sin(x)', 'y*x']

print(p.form)
print(p.rhs)
print(p.f)
print(p.uB)
print(p.uE)
print(p.wind)

lvp = LinearVariationalProblem(p.form, p.rhs, p.sol, p.bcs)
lvs = LinearVariationalSolver(lvp, solver_parameters={'ksp_type': 'cg',
                                                      'pc_type': 'hypre'})
lvs.solve()
print (np.linalg.norm(p.sol.dat.data[:, 0] - Expression(p.uE, V).dat.data[:, 0]))
