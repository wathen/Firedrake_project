from .BaseClass import BaseClass
from .Poisson import Poisson
import firedrake as fd
import sympy


class ConvecDiff(Poisson):

    def __init__(self, V, nu, wind=0, opt='linear'):
        super(ConvecDiff, self).__init__(V, nu)
        self.nu = nu
        self._opt = opt
        self.wind_ = wind
        if opt == 'linear':
            self._wind = self.Expression(wind)
            self.wind_ = wind
        else:
            self._wind = super().sol
            self.wind_ = uE

    def advection(self):
        if self._opt == 'linear':
            return fd.inner((fd.grad(self.sol) * self._wind), self.v) * fd.dx
        elif self._opt == 'nonlinear':
            return fd.inner((fd.grad(self._wind) * self._wind), self.v) * fd.dx

    @property
    def form(self):
        if not self._form:
            self._form = self.laplacian() + self.advection()
        return self._form

    def Advetion_exact_differentiator(self, u):
        element_dim = len(u)
        while len(u) < 3:
            u.append('0')
            self.wind_.append('0')

        x, y, z = sympy.symbols('x y z')
        u = [item.replace('math', 'sympy') for item in u]

        if element_dim >= 1:
            L1 = eval(self.wind_[0]) * sympy.diff(eval(u[0]), x) + eval(self.wind_[1]) * sympy.diff(eval(u[0]), y) + eval(self.wind_[1]) * sympy.diff(eval(u[0]), z)
            out = [str(L1)]
        if element_dim >= 2:
            L2 = eval(self.wind_[0]) * sympy.diff(eval(u[1]), x) + eval(self.wind_[1]) * sympy.diff(eval(u[1]), y) + eval(self.wind_[1]) * sympy.diff(eval(u[1]), z)
            out.append(str(L2))
        if element_dim >= 3:
            L3 = eval(self.wind_[0]) * sympy.diff(eval(u[2]), x) + eval(self.wind_[1]) * sympy.diff(eval(u[2]), y) + eval(self.wind_[1]) * sympy.diff(eval(u[2]), z)
            out.append(str(L3))

        return super().convert_sympy_firedrake_expresion(out)

    @property
    def bcs(self):
        if not self._bcs:
            if self._uE:
                self._uB = self._uE
            self._bcs = [fd.DirichletBC(self._V,
                                        self.Expression(self._uB),
                                        'on_boundary')]
        return self._bcs

    @property
    def rhs(self):
        if not self._f:
            if self._uE:
                self._f = self.laplacian_exact_differentiator(self._uE) \
                    + self.Advetion_exact_differentiator(self._uE)
            self._rhs = fd.inner(self.Expression(self._f), self.v) * fd.dx
        return self._rhs

    @property
    def wind(self):
        return self._wind

    @wind.setter
    def wind(self, value):
        self._wind = value
