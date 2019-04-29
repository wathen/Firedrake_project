from .Poisson import Poisson
from firedrake import inner, dx, dot, grad, DirichletBC, FacetNormal, nabla_grad, VectorFunctionSpace


class Helmholtz(Poisson):

    def __init__(self, V, nu, opt={'form': 'linear'}):
        super(Helmholtz, self).__init__(V, nu)

    def mass(self):
        if self.opt['form'] == 'linear' or self.opt['form'] == 'nonlinear':
            return self.k * inner(self.sol, self.v) * dx
        elif self.opt['form'] == 'bilinear':
            return self.k * inner(self.u, self.v) * dx

    @property
    def form(self):
        if not self._form:
            self._form = -self.laplacian() + self.mass()
        return self._form

    @property
    def f(self):
        if not self._f:
            if self.uE:
                laplacian = self.ExactSol.laplacian(self.uE)
                self._f = self.ExactSol.vec_add(laplacian, self.uE, self.k)
        return self._f

    @property
    def k(self):
        return self._k

    @k.setter
    def k(self, value):
        self._k = value
