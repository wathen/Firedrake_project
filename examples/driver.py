from firedrake import *
from project import ConvecDiff
from project.utils import Expression


class Lapl_example(ConvecDiff):

    def __init__(self, V, nu):
        super(Lapl_example, self).__init__(V, nu)


mesh = UnitCubeMesh(2, 2, 2)
V = VectorFunctionSpace(mesh, 'CG', 1)
p = ConvecDiff(V, 1, ['y**3', 'y**3'])
p.uE = ['y**3', 'y**3', 'y**3']


print(p.bcs)
print(p.rhs)
print(p.form)
print(p.sol)
#
lvp = LinearVariationalProblem(p.form, p.rhs, p.sol, p.bcs)
lvs = LinearVariationalSolver(lvp)

lvs.solve()
print (p.sol.dat.data)
