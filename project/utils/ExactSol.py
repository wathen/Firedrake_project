from sympy import *
import sympy.vector as spv
import re


class ExactSol(object):

    def vector_gradient(self, expr):
        R = spv.CoordSys3D('R')
        x, y, z = R.base_scalars()

        grad = []
        for i in range(len(expr)):
            grad.append(spv.gradient(eval(expr[i])))
        return self.post_process(grad)

    def laplacian(self, expr):
        R = spv.CoordSys3D('R')
        x, y, z = R.base_scalars()
        out = []
        for i in range(len(expr)):
            out.append(spv.divergence(spv.gradient(eval(expr[i]))))
        return self.post_process(out)

    def advection(self, expr, wind):
        R = spv.CoordSys3D('R')
        x, y, z = R.base_scalars()
        grad = []
        for i in range(len(expr)):
            grad.append(spv.gradient(eval(expr[i])))
        out = []
        for g in grad:
            g_list = self.ijk_to_list(g)
            temp = ''
            for (w, g) in zip(wind, g_list):
                if w == '0' or g == '0':
                    temp += '+0'
                else:
                    temp += '+(' + w + ')*(' + g + ')'
            out.append(temp)

        return out

    def curl(self, expr):
        R = spv.CoordSys3D('R')
        x, y, z = R.base_scalars()
        expr_len = len(expr)
        if expr_len == 1:
            expr = ['0'] + ['0'] + expr
        elif expr_len == 2:
            expr = expr + ['0']

        expr_ijk = self.list_to_ijk(expr)
        curl_expr = spv.curl(eval(expr_ijk))
        curl_expr_list = self.ijk_to_list(curl_expr)
        if expr_len == 1:
            a, b, *_ = curl_expr_list
            return [a, b]
        elif expr_len == 2:
            *_, a = curl_expr_list
            return [a]
        else:
            return curl_expr_list

    def curlcurl(self, expr):
        return self.curl(self.curl(expr))

    def cross(self, a, b):
        R = spv.CoordSys3D('R')
        x, y, z = R.base_scalars()
        a_len = len(a)
        if a_len == 1:
            a = ['0'] + ['0'] + a
        elif a_len == 2:
            a = a + ['0']

        b_len = len(b)
        if b_len == 1:
            b = ['0'] + ['0'] + b
        elif b_len == 2:
            b = b + ['0']

        a_ijk = self.list_to_ijk(a)
        b_ijk = self.list_to_ijk(b)
        cross_expr = spv.cross(eval(a_ijk), eval(b_ijk))
        cross_expr_list = self.ijk_to_list(cross_expr)
        return cross_expr_list

    def post_process(self, a):
        out = str(a).replace('R.', '')
        out = out.replace('[', '').replace(']', '').split(',')
        return out

    def ijk_to_list(self, text, pattern=['R.i', 'R.j', 'R.k']):
        R = spv.CoordSys3D('R')
        x, y, z = R.base_scalars()
        r = eval(str(text))
        if isinstance(r, int):
            l = [r, r, r]
        else:
            l = [r.coeff(R.i), r.coeff(R.j), r.coeff(R.k)]
        l = [str(x).replace('R.', '') for x in l]
        return l

    def list_to_ijk(self, coord):
        out = ''
        y = 0
        c = ['i', 'j', 'k']
        for i in coord:
            if i != '0':
                out += " + (" + str(i) + ")*R." + c[y]
            y += 1
        return out

    def scalar_mult(self, scalar, expr):
        out = []
        for i in expr:
            out.append(str(scalar) + '*' + i)
        return out
