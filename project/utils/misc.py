from firedrake import *


def Expression(f, V):
    mesh_dim = V.mesh().cell_dimension()
    f_len = len(f)
    if mesh_dim == 1:
        x = SpatialCoordinate(V.mesh())
    elif mesh_dim == 2:
        (x, y) = SpatialCoordinate(V.mesh())
    elif mesh_dim == 3:
        (x, y, z) = SpatialCoordinate(V.mesh())

    if isinstance(f, list):
        if f_len == 1:
            out = interpolate(eval(f[0]), V)
        elif f_len == 2:
            out = interpolate(as_vector([eval(f[0]), eval(f[1])]), V)
        elif f_len == 3:
            out = interpolate(as_vector([eval(f[0]), eval(f[1]), eval(f[2])]), V)
        else:
            raise RuntimeError('Input list for function value is of '
                               'dimension {}, need to be at most a 3D '
                               'vector'.format(len(f)))
        return out
    else:
        raise RuntimeError('Input function f needs to be list')
