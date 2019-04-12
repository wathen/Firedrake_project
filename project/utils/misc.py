from firedrake import Function, SpatialCoordinate, as_vector, interpolate


def Expression(f, mesh):
    x = SpatialCoordinate(mesh)
    if isinstance(f, list):
        if len(f) == 1:
            out = eval(f[0])
        elif len(f) == 2:
            out = as_vector([eval(f[0]), eval(f[1])])
        elif len(f) == 3:
            out = as_vector([eval(f[0]), eval(f[1]), eval(f[3])])
        else:
            raise RuntimeError('Input list for function value is of '
                               'dimension {}, need to be at most a 3D '
                               'vector'.format(len(f)))
        return out
    else:
        raise RuntimeError('Input function f needs to be list')
