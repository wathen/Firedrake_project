from firedrake import FunctionSpace, VectorFunctionSpace, SpatialCoordinate, interpolate, as_vector, Constant, exp, sin, cos, tan, asin, acos, atan, atan_2, cosh, sinh, tanh, pi


def Expression(f, V, degree_raise=None):
    if degree_raise != None:
        V = FunctionSpaceDegreeRaise(V, degree_raise)

    mesh_dim = V.mesh().cell_dimension()
    f_len = len(f)

    if mesh_dim == 1:
        x = SpatialCoordinate(V.mesh())
    elif mesh_dim == 2:
        (x, y) = SpatialCoordinate(V.mesh())
    elif mesh_dim == 3:
        (x, y, z) = SpatialCoordinate(V.mesh())

    f_eval = [None] * f_len
    for i in range(f_len):
        if 'x' not in f[i] and 'y' not in f[i] and 'z' not in f[i]:
            f_eval[i] = Constant(eval(f[i]))
        else:
            f_eval[i] = eval(f[i])

    if isinstance(f, list):
        if f_len == 1:
            out = interpolate(f_eval[0], V)
        elif f_len == 2:
            out = interpolate(as_vector([f_eval[0], f_eval[1]]), V)
        elif f_len == 3:
            out = interpolate(as_vector([f_eval[0], f_eval[1], f_eval[2]]), V)
        else:
            raise RuntimeError('Input list for function value is of '
                               'dimension {}, need to be at most a 3D '
                               'vector'.format(len(f)))
        return out
    else:
        raise RuntimeError('Input function f needs to be list')


def FunctionSpaceDegreeRaise(V, degree_raise):
    element = V._ufl_element

    degree = element.degree(0)
    space_type = repr(element)[0]
    element_type = element.family()

    mesh = V.mesh()
    if space_type == 'V':
        V = VectorFunctionSpace(mesh, element_type, degree + degree_raise)
    elif space_type == 'F':
        V = FunctionSpace(mesh, element_type, degree + degree_raise)

    return V
