import os
import sys
import unittest
parent_dir = os.path.dirname(os.getcwd())
sys.path.append(parent_dir)

import ExactSol

exact_sol = ExactSol.ExactSol()


class ExactSolTests_ijk_to_list(unittest.TestCase):

    def test_ijk_to_list_1d(self):
        expr = '3*z**2*R.k'
        expected = ['0', '0', '3*z**2']
        assert expected == exact_sol.ijk_to_list(expr)

    def test_ijk_to_list_2d(self):
        expr = 'R.i + 3*z**2*R.k'
        expected = ['1', '0', '3*z**2']
        assert expected == exact_sol.ijk_to_list(expr)

    def test_ijk_to_list_3d(self):
        expr = 'R.j + R.i + 3*z**2*R.k'
        expected = ['1', '1', '3*z**2']
        assert expected == exact_sol.ijk_to_list(expr)


class ExactSolTests_list_to_ijk(unittest.TestCase):

    def test_list_to_ijk_1d(self):
        expr = ['tanh(x+5*z)']
        expected = ' + (tanh(x+5*z))*R.i'
        assert expected == exact_sol.list_to_ijk(expr)

    def test_list_to_ijk_2d(self):
        expr = ['5', 'sin(x+y)']
        expected = ' + (5)*R.i + (sin(x+y))*R.j'
        assert expected == exact_sol.list_to_ijk(expr)

    def test_list_to_ijk_3d(self):
        expr = ['0', '0', '3*z**2']
        expected = ' + (3*z**2)*R.k'
        assert expected == exact_sol.list_to_ijk(expr)


class ExactSolTests_laplacian(unittest.TestCase):

    def test_laplacian_1d(self):
        expr = ['x**2+sin(y)']
        expectedLaplacian = ['-sin(y) + 2']
        assert expectedLaplacian == exact_sol.laplacian(expr)

    def test_laplacian_2d(self):
        expr = ['x**2+sin(y)', 'exp(y)+cos(z)']
        expectedLaplacian = ['-sin(y) + 2', ' exp(y) - cos(z)']
        assert expectedLaplacian == exact_sol.laplacian(expr)

    def test_laplacian_3d(self):
        expr = ['x**2+sin(y)', 'exp(y)+cos(z)', 'tan(z+x)']
        expectedLaplacian = ['-sin(y) + 2',
                             ' exp(y) - cos(z)',
                             ' 2*(2*tan(x + z)**2 + 2)*tan(x + z)']
        assert expectedLaplacian == exact_sol.laplacian(expr)


class ExactSolTests_vector_gradient(unittest.TestCase):

    def test_vector_gradient_1d(self):
        expr = ['x**2+sin(y)']
        expectedvector_gradient = ['2*x*i + (cos(y))*j']
        assert expectedvector_gradient == exact_sol.vector_gradient(expr)

    def test_vector_gradient_2d(self):
        expr = ['x**2+sin(y)', 'exp(y)+cos(z)']
        expectedvector_gradient = ['2*x*i + (cos(y))*j',
                                   ' (exp(y))*j + (-sin(z))*k']
        assert expectedvector_gradient == exact_sol.vector_gradient(expr)

    def test_vector_gradient_3d(self):
        expr = ['x**2+sin(y)', 'exp(y)+cos(z)', 'tan(z+x)']
        expectedvector_gradient = ['2*x*i + (cos(y))*j',
                                   ' (exp(y))*j + (-sin(z))*k',
                                   ' (tan(x + z)**2 + 1)*i + (tan(x + z)**2 + 1)*k']
        assert expectedvector_gradient == exact_sol.vector_gradient(expr)


class ExactSolTests_advection(unittest.TestCase):

    def test_advection_1d(self):
        expr = ['x**2']
        wind = ['tan(y)']
        expectedadvection = ['+(tan(y))*(2*x)']
        assert expectedadvection == exact_sol.advection(expr, wind)

    def test_advection_2d(self):
        expr = ['x**2+sin(y)', 'exp(y)+cos(z)']
        wind = ['tan(y)', 'exp(x)']
        expectedadvection = ['+(tan(y))*(2*x)+(exp(x))*(cos(y))',
                             '+(exp(x))*(exp(y))']
        assert expectedadvection == exact_sol.advection(expr, wind)

    def test_advection_3d(self):
        expr = ['x**2+sin(y)', 'exp(y)+cos(z)', 'tan(z+x)']
        wind = ['tan(y)', 'exp(x)', 'x*exp(z)']
        expectedadvection = ['+(tan(y))*(2*x)+(exp(x))*(cos(y))',
                             '+(exp(x))*(exp(y))+(x*exp(z))*(-sin(z))',
                             '+(tan(y))*(tan(x + z)**2 + 1)+(x*exp(z))*(tan(x + z)**2 + 1)']
        assert expectedadvection == exact_sol.advection(expr, wind)


class ExactSolTests_curl(unittest.TestCase):

    def test_curl_1d(self):
        expr = ['x**2+sin(y)']
        expectedcurl = ['cos(y)', '-2*x']
        assert expectedcurl == exact_sol.curl(expr)

    def test_curl_2d(self):
        expr = ['x**2+sin(y)', 'exp(y)+cos(z)']
        expectedcurl = ['-cos(y)']
        assert expectedcurl == exact_sol.curl(expr)

    def test_curl_3d(self):
        expr = ['x**2+sin(y)', 'exp(y)+cos(z)', 'tan(z+x)']
        expectedcurl = ['sin(z)', '-tan(x + z)**2 - 1', '-cos(y)']
        assert expectedcurl == exact_sol.curl(expr)
