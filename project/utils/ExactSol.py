from sympy import *
import sympy.vector as spv
import re

R = spv.CoordSys3D('R')
x, y, z = R.base_scalars()


class ExactSol(object):

    def vector_gradient(self, expr):

        grad = []
        for i in range(len(expr)):
            print(self.ijk_to_list(spv.gradient(eval(expr[i]))))
            grad.append(self.ijk_to_list(spv.gradient(eval(expr[i]))))
        return grad

    def laplacian(self, expr):
        out = []
        for i in range(len(expr)):
            out.append(spv.divergence(spv.gradient(eval(expr[i]))))
        return self.post_process(out)

    def advection(self, expr, wind):
        grad = self.vector_gradient(expr)
        out = []
        for g_list in grad:
            temp = ''
            for (w, g) in zip(wind, g_list):
                if w == '0' or g == '0':
                    temp += '+0'
                else:
                    temp += '+(' + w + ')*(' + g + ')'
            out.append(temp)

        return out

    def curl(self, expr):
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
            if i == '0':
                out.append(i)
            else:
                out.append(str(scalar) + '*' + i)
        return out

    def vec_add(self, expr1, expr2, alpha=1):
        out = []
        for (a, b) in zip(expr1, expr2):
            if a == '0' and b == '0':
                out.append('0')
            else:
                out.append(a + '+' + str(alpha) + '*(' + b + ')')

        return out
