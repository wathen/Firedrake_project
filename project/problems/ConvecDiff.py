from .BaseClass import BaseClass
from .Poisson import Poisson
from firedrake import inner, dx, dot, grad, DirichletBC, FacetNormal, nabla_grad, VectorFunctionSpace
from project.utils import Expression


class ConvecDiff(Poisson):

    def __init__(self, V, nu, opt={'form': 'linear'}):
        super(ConvecDiff, self).__init__(V, nu)

    def advection(self):
        n = FacetNormal(self.mesh)
        if len(self.wind) > 1:
            element = self.V._ufl_element

            degree = element.degree(0)
            space_type = repr(element)[0]
            element_type = element.family()

            mesh = self.V.mesh()
            V = VectorFunctionSpace(mesh, element_type, degree)
        else:
            V = self.V
        if self.opt['form'] == 'linear':
            w_ = Expression(self.wind, V)
            u_ = self.sol
        elif self.opt['form'] == 'bilinear':
            w_ = Expression(self.wind, V)
            u_ = self.u
        elif self.opt['form'] == 'nonlinear':
            w_ = self.sol
            u_ = self.sol
        return inner(dot(w_, nabla_grad(u_)), self.v) * dx

    @property
    def form(self):
        if not self._form:
            self._form = self.laplacian() + self.advection()
        return self._form

    @property
    def f(self):
        if not self._f:
            if self.uE:
                laplacian = self.ExactSol.laplacian(self.uE)
                advection = self.ExactSol.advection(self.uE, self.wind)
                laplacian = self.ExactSol.scalar_mult(-self.nu, laplacian)
                self._f = [x[0] + ' + ' + x[1]
                           for x in zip(laplacian, advection)]
        return self._f

    @property
    def wind(self):
        return self._wind

    @wind.setter
    def wind(self, value):
        self._wind = value
