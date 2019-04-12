from .BaseClass import BaseClass
import firedrake as fd
import sympy


class Poisson(BaseClass):

    def __init__(self, V, nu):
        super(Poisson, self).__init__(V)
        self.nu = nu

    def laplacian(self):
        return self.nu * fd.inner(fd.grad(self.sol), fd.grad(self.v)) * fd.dx

    @property
    def form(self):
        if not self._form:
            self._form = self.laplacian()
        return self._form

    def laplacian_exact_differentiator(self, u):
        element_dim = len(u)
        x, y, z = sympy.symbols('x y z')
        u = [item.replace('math', 'sympy') for item in u]
        if element_dim >= 1:
            L1 = sympy.diff(eval(u[0]), x, x) \
                + sympy.diff(eval(u[0]), y, y) \
                + sympy.diff(eval(u[0]), z, z)
            out = [str(L1)]
        if element_dim >= 2:
            L2 = sympy.diff(eval(u[1]), x, x) \
                + sympy.diff(eval(u[1]), y, y) \
                + sympy.diff(eval(u[1]), z, z)
            out.append(str(L2))
        if element_dim >= 3:
            L3 = sympy.diff(eval(u[2]), x, x) \
                + sympy.diff(eval(u[2]), y, y) \
                + sympy.diff(eval(u[2]), z, z)
            out.append(str(L3))
        out = [str(self.nu) + ' * (' + str(item) + ')' for item in out]
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
                self._f = self.laplacian_exact_differentiator(self._uE)
            print (self.Expression(self._f))
            self._rhs = fd.inner(self.Expression(self._f), self.v) * fd.dx(self.mesh)
        return self._rhs
