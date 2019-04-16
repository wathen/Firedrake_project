from .BaseClass import BaseClass
import firedrake as fd
from sympy import *
from project.utils import Expression


class Poisson(BaseClass):

    def __init__(self, V, nu, opt={'form': 'linear'}):
        super(Poisson, self).__init__(V)
        self.nu = nu
        self.opt = opt

    def laplacian(self):
        if self.opt['form'] == 'linear':
            return self.nu * fd.inner(fd.grad(self.sol), fd.grad(self.v)) * fd.dx
        elif self.opt['form'] == 'bilinear':
            return self.nu * fd.inner(fd.grad(self.u), fd.grad(self.v)) * fd.dx

    @property
    def form(self):
        if not self._form:
            self._form = self.laplacian()
        return self._form

    def laplacian_exact_differentiator(self, u):
        element_dim = len(u)
        x, y, z = symbols('x y z')
        u = [item.replace('math', 'sympy') for item in u]
        if element_dim >= 1:
            L1 = diff(eval(u[0]), x, x) \
                + diff(eval(u[0]), y, y) \
                + diff(eval(u[0]), z, z)
            out = [str(-L1)]
        if element_dim >= 2:
            L2 = diff(eval(u[1]), x, x) \
                + diff(eval(u[1]), y, y) \
                + diff(eval(u[1]), z, z)
            out.append(str(-L2))
        if element_dim >= 3:
            L3 = diff(eval(u[2]), x, x) \
                + diff(eval(u[2]), y, y) \
                + diff(eval(u[2]), z, z)
            out.append(str(-L3))
        out = [str(self.nu) + ' * (' + str(item) + ')' for item in out]
        return super().convert_sympy_firedrake_expresion(out)

    @property
    def bcs(self):
        if not self._bcs:
            self._bcs = [fd.DirichletBC(self._V,
                                        Expression(self.uB, self.V),
                                        'on_boundary')]
        return self._bcs

    @property
    def f(self):
        if not self._f:
            if self._uE:
                self._f = self.laplacian_exact_differentiator(self.uE)
        return self._f

    @property
    def rhs(self):
        if not self._rhs:
            self._rhs = fd.inner(Expression(self.f, self.V), self.v) * fd.dx
        return self._rhs
