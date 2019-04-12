from firedrake import FunctionSpace, VectorFunctionSpace


def TaylorHood(mesh, order):

    V = VectorFunctionSpace(mesh, 'CG', order + 1)
    Q = FunctionSpace(mesh, 'CG', order)

    return V * Q, {'V': V, 'P': Q}


def MixedNedelec(mesh, order):

    V = VectorFunctionSpace(mesh, 'N1curl', order)
    Q = FunctionSpace(mesh, 'CG', order)

    return V * Q, {'M': V, 'L': Q}


def TaylorHoodMixedNedelec(mesh, order):

    TH, TH_spaces = TaylorHood(mesh, order)
    MN, MN_spaces = MixedNedelec(mesh, order)

    return TH * MN, dict(TH_spaces, **MN_spaces)
