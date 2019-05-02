from .BaseClass import BaseClass
from firedrake import inner, dx, grad, div, Expression, DirichletBC
from project.utils import Expression


class Stokes(BaseClass):
    """
    Definition of the steady state Poisson class which implements the finite
    element solution for Poissons equation:
    .. math::
        -\\Delta u + \\nabla p = f on \\Omega
        -\\nabla \\cdot u      = 0 on \\Omega
        with u = u_b on \\partial\\Omega
    """

    def __init__(self, V, nu, opt={'form': 'linear'}):
        V_mixed, V_dict = V
        super(Stokes, self).__init__(V_mixed)
        self.V_dict = V_dict
        self.nu = nu
        self.opt = opt
        self._pE = None

    def stokes(self):
        v_ = self.v
        if self.opt['form'] == 'linear' or self.opt['form'] == 'nonlinear':
            u_ = self.sol.split()
        elif self.opt['form'] == 'bilinear':
            u_ = self.u
        l = self.nu * inner(grad(u_[0]), grad(v_[0])) * dx
        b = -div(v_[0]) * u_[1] * dx
        bt = -div(u_[0]) * v_[1] * dx
        return l + b + bt

    @property
    def form(self):
        if not self._form:
            self._form = self.stokes()
        return self._form

    @property
    def f(self):
        if not self._f:
            if self._uE:
                laplacian = self.ExactSol.laplacian(self.uE)
                grad_p = self.ExactSol.vector_gradient(self.pE)[0]
                self._f = self.ExactSol.vec_add(grad_p, laplacian, -self.nu)
        return self._f

    @property
    def rhs(self):
        if not self._rhs:
            f_ = Expression(self.f, self.V_dict['V'])
            print(f_)
            self._rhs = inner(f_, self.v[0]) * dx

        return self._rhs

    @property
    def bcs(self):
        if not self._bcs:
            self._bcs = [DirichletBC(self.V.sub(0),
                                     Expression(self.uB, self.V_dict['V']),
                                     'on_boundary')]
        return self._bcs

    @property
    def pE(self):
        return self._pE

    @pE.setter
    def pE(self, value):
        self._pE = value
