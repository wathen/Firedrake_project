from .BaseClass import BaseClass
from .Poisson import Poisson
import firedrake as fd
from sympy import *
from project.utils import Expression


class ConvecDiff(Poisson):

    def __init__(self, V, nu, opt={'form': 'linear'}):
        super(ConvecDiff, self).__init__(V, nu)

    def advection(self):
        w = Expression(self.wind, self.V)
        if self.opt['form'] == 'linear':
            return fd.inner((fd.grad(self.u) * w), self.v) * fd.dx
        if self.opt['form'] == 'bilinear':
            return fd.inner((fd.grad(w) * w), self.v) * fd.dx

    @property
    def form(self):
        if not self._form:
            self._form = self.laplacian() + self.advection()
        return self._form

    def advection_exact_differentiator(self, u, w):
        element_dim = len(u)
        w_len = len(w)
        # while len(w) < 3:
        #     if w_len == 1:
        #         w.append(w[0])
        #     else:
        #         w.append('0')
        mesh_dim = self.V.mesh().cell_dimension()
        if mesh_dim == 1:
            x = symbols('x')
        elif mesh_dim == 2:
            x, y = symbols('x y')
        elif mesh_dim == 3:
            x, y, z = symbols('x y z')
        print([x, y])

        if element_dim >= 1:
            A1 = eval(w[0]) * diff(eval(u[0]), x) + eval(w[1]) * diff(eval(u[0]), y) + eval(w[1]) * diff(eval(u[0]), z)
            out = [str(A1)]
        if element_dim >= 2:
            A2 = eval(w[0]) * diff(eval(u[1]), x) + eval(w[1]) * diff(eval(u[1]), y) + eval(w[1]) * diff(eval(u[1]), z)
            out.append(str(A2))
        if element_dim >= 3:
            A3 = eval(w[0]) * diff(eval(u[2]), x) + eval(w[1]) * diff(eval(u[2]), y) + eval(w[1]) * diff(eval(u[2]), z)
            out.append(str(A3))
        return super().convert_sympy_firedrake_expresion(out)

    def convect(self, w, u, x):
        return eval(w) * diff(eval(u), x)

    @property
    def f(self):
        if not self._f:
            if self.uE:
                L = self.laplacian_exact_differentiator(self.uE)
                A = self.advection_exact_differentiator(self.uE, self.wind)
                self._f = [x[0] + ' + ' + x[1] for x in zip(L, A)]
        return self._f

    @property
    def rhs(self):
        if not self._rhs:
            self._rhs = fd.inner(Expression(self.f, self.V), self.v) * fd.dx
        return self._rhs

    @property
    def bcs(self):
        if not self._bcs:
            self._bcs = [fd.DirichletBC(self.V,
                                        Expression(self.uB, self.V),
                                        'on_boundary')]
        return self._bcs

    @property
    def wind(self):
        return self._wind

    @wind.setter
    def wind(self, value):
        self._wind = value
