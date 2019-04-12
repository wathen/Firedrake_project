from firedrake import TestFunction, TrialFunction, \
    TestFunctions, TrialFunctions, SpatialCoordinate, \
    as_vector, Function, interpolate
from abc import ABCMeta, abstractmethod
import math


class BaseClass(object):

    __metaclass__ = ABCMeta

    def __init__(self, V):
        self._V = V
        self.mesh = self._V.mesh()
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

    def Expression(self, f):
        mesh_dim = self.mesh.cell_dimension()

        if mesh_dim == 1:
            x = SpatialCoordinate(self._V.mesh())
        elif mesh_dim == 2:
            (x, y) = SpatialCoordinate(self._V.mesh())
        elif mesh_dim == 3:
            (x, y, z) = SpatialCoordinate(self._V.mesh())

        if isinstance(f, list):
            if len(f) == 1:
                out = interpolate(eval(f[0]), self._V)
            elif len(f) == 2:
                out = interpolate(as_vector([eval(f[0]), eval(f[1])]), self._V)
            elif len(f) == 3:
                out = interpolate(as_vector([eval(f[0]), eval(f[1]), eval(f[2])]), self._V)
            else:
                raise RuntimeError('Input list for function value is of '
                                   'dimension {}, need to be at most a 3D '
                                   'vector'.format(len(f)))
            return out
        else:
            raise RuntimeError('Input function f needs to be list')

    def convert_sympy_firedrake_expresion(self, u):
        u = [item.replace('pi', 'math.pi') for item in u]
        u = [item.replace('cos', 'math.cos') for item in u]
        u = [item.replace('sin', 'math.sin') for item in u]
        u = [item.replace('tan', 'math.tan') for item in u]
        u = [item.replace('exp', 'math.exp') for item in u]
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
