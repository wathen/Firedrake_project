from .BaseClass import BaseClass
from firedrake import inner, dx, grad


class Poisson(BaseClass):
    """
    Definition of the steady state Poisson class which implements the finite
    element solution for Poissons equation:
    .. math::
        -\\Delta u = f on \\Omega
        with u = u_b on \\partial\\Omega

    Types of data:
        - strings: f, uE, uB
        - firedake operators: form, rhs, bcs, sol
    """

    def __init__(self, V, nu, opt={'form': 'linear'}):
        super(Poisson, self).__init__(V)
        self.nu = nu
        self.opt = opt

    def laplacian(self):
        if self.opt['form'] == 'linear' or self.opt['form'] == 'nonlinear':
            return self.nu * inner(grad(self.sol), grad(self.v)) * dx
        elif self.opt['form'] == 'bilinear':
            return self.nu * inner(grad(self.u), grad(self.v)) * dx

    @property
    def form(self):
        if not self._form:
            self._form = self.laplacian()
        return self._form

    @property
    def f(self):
        if not self._f:
            if self._uE:
                laplacian = self.ExactSol.laplacian(self.uE)
                self._f = self.ExactSol.scalar_mult(-self.nu, laplacian)
        return self._f
