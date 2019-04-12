from firedrake import *
from project import Poisson
from project import ConvecDiff
from project.utils import Expression


n = 2**2.
mesh = UnitSquareMesh(n, n)
V = VectorFunctionSpace(mesh, 'CG', 1)

# u = TestFunction(V)
# v = TrialFunction(V)

# (x, y) = SpatialCoordinate(mesh)
# f = interpolate(as_vector([0 * x, 0 * x]), V)
# print (f)
# inner(f, v) * dx

p = ConvecDiff(V, 1, wind=['x*x', 'y*x'], opt='linear')
# p.wind = ['x*x', 'y*x']
p.uE = ['x', 'y']

print(p.bcs)
# print(p.rhs)
print(p.form - p.rhs)
print(p.sol)
#
lvp = NonlinearVariationalProblem(p.form - p.rhs, p.sol, p.bcs)
lvs = NonlinearVariationalSolver(lvp)

lvs.solve()
print (p.sol.dat.data)
print (p.wind.dat.data)
print (p.Expression(p.uE).dat.data)
print(p.uE)
