from firedrake import *
from abc import ABCMeta, abstractmethod
import math


class BaseClass(object):

    __metaclass__ = ABCMeta

    def __init__(self, V):
        self._V = V
        self.mesh = self.V.mesh()
        self.test_trial_functions(self._V)
        self._form = None
        self._rhs = None
        self._bcs = None
        self._f = None
        self._uB = None
        self._uE = None
        self._sol = None

    @property
    def dim(self):
        return self._V.dim()

    @property
    def V(self):
        return self._V

    @V.setter
    def V(self, value):
        self._V = value

    @property
    def sol(self):
        if not self._sol:
            self._sol = Function(self._V)
        return self._sol

    @sol.setter
    def sol(self, value):
        self._sol = value

    @property
    def form(self):
        return self._form

    @form.setter
    def form(self, value):
        self._form = value

    @property
    def rhs(self):
        return self._rhs

    @rhs.setter
    def rhs(self, value):
        self._rhs = value

    @property
    def bcs(self):
        return self._bcs

    @bcs.setter
    def bcs(self, value):
        self._bcs = value

    @property
    def f(self):
        return self._f

    @f.setter
    def f(self, value):
        self._f = value

    @property
    def uB(self):
        if not self._uB:
            if self._uE:
                self._uB = self.uE
        return self._uB

    @uB.setter
    def uB(self, value):
        self._uB = value

    @property
    def uE(self):
        return self._uE

    @uE.setter
    def uE(self, value):
        self._uE = value

    def test_trial_functions(self, V):
        try:
            dim = W.num_sub_spaces()
        except:
            dim = 1

        if dim == 1:
            self.u = TrialFunction(V)
            self.v = TestFunction(V)
        else:
            self.u = TrialFunctions(V)
            self.v = TestFunctions(V)

    def convert_sympy_firedrake_expresion(self, u):
        u = [item.replace('pi', 'math.pi') for item in u]
        return u

    def has_nullspace(self):
        raise NotImplementedError

    def nullspace(self, V):
        if self.has_nullspace():
            MVSB = MixedVectorSpaceBasis
            return MVSB(V, [V.sub(0), VectorSpaceBasis(constant=True)])
        else:
            return None

    def mesh_size(self, u):
        return CellSize(u.ufl_domain())
