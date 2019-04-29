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
        expect_ = ['0', '0', '3*z**2']
        assert expect_ == exact_sol.ijk_to_list(expr)

    def test_ijk_to_list_2d(self):
        expr = 'R.i + 3*z**2*R.k'
        expect_ = ['1', '0', '3*z**2']
        assert expect_ == exact_sol.ijk_to_list(expr)

    def test_ijk_to_list_3d(self):
        expr = 'R.j + R.i + 3*z**2*R.k'
        expect_ = ['1', '1', '3*z**2']
        assert expect_ == exact_sol.ijk_to_list(expr)


class ExactSolTests_list_to_ijk(unittest.TestCase):

    def test_list_to_ijk_1d(self):
        expr = ['tanh(x+5*z)']
        expect_ = ' + (tanh(x+5*z))*R.i'
        assert expect_ == exact_sol.list_to_ijk(expr)

    def test_list_to_ijk_2d(self):
        expr = ['5', 'sin(x+y)']
        expect_ = ' + (5)*R.i + (sin(x+y))*R.j'
        assert expect_ == exact_sol.list_to_ijk(expr)

    def test_list_to_ijk_3d(self):
        expr = ['0', '0', '3*z**2']
        expect_ = ' + (3*z**2)*R.k'
        assert expect_ == exact_sol.list_to_ijk(expr)


class ExactSolTests_laplacian(unittest.TestCase):

    def test_laplacian_1d(self):
        expr = ['x**2+sin(y)']
        expect_Laplacian = ['-sin(y) + 2']
        assert expect_Laplacian == exact_sol.laplacian(expr)

    def test_laplacian_2d(self):
        expr = ['x**2+sin(y)', 'exp(y)+cos(z)']
        expect_Laplacian = ['-sin(y) + 2', ' exp(y) - cos(z)']
        assert expect_Laplacian == exact_sol.laplacian(expr)

    def test_laplacian_3d(self):
        expr = ['x**2+sin(y)', 'exp(y)+cos(z)', 'tan(z+x)']
        expect_Laplacian = ['-sin(y) + 2',
                            ' exp(y) - cos(z)',
                            ' 2*(2*tan(x + z)**2 + 2)*tan(x + z)']
        assert expect_Laplacian == exact_sol.laplacian(expr)


class ExactSolTests_vector_gradient(unittest.TestCase):

    def test_vector_gradient_1d(self):
        expr = ['x**2+sin(y)']
        expect_vector_gradient = [['2*x', 'cos(y)', '0']]
        assert expect_vector_gradient == exact_sol.vector_gradient(expr)

    def test_vector_gradient_2d(self):
        expr = ['x**2+sin(y)', 'exp(y)+cos(z)']
        expect_vector_gradient = [['2*x', 'cos(y)', '0'],
                                  ['0', 'exp(y)', '-sin(z)']]

        assert expect_vector_gradient == exact_sol.vector_gradient(expr)

    def test_vector_gradient_3d(self):
        expr = ['x**2+sin(y)', 'exp(y)+cos(z)', 'tan(z+x)']
        expect_vector_gradient = [['2*x', 'cos(y)', '0'],
                                  ['0', 'exp(y)', '-sin(z)'],
                                  ['tan(x + z)**2 + 1', '0', 'tan(x + z)**2 + 1']]

        assert expect_vector_gradient == exact_sol.vector_gradient(expr)


class ExactSolTests_advection(unittest.TestCase):

    def test_advection_1d(self):
        expr = ['x**2']
        wind = ['tan(y)']
        expect_advection = ['+(tan(y))*(2*x)']
        assert expect_advection == exact_sol.advection(expr, wind)

    def test_advection_2d(self):
        expr = ['x**2+sin(y)', 'exp(y)+cos(z)']
        wind = ['tan(y)', 'exp(x)']
        expect_advection = ['+(tan(y))*(2*x)+(exp(x))*(cos(y))',
                            '+0+(exp(x))*(exp(y))']
        assert expect_advection == exact_sol.advection(expr, wind)

    def test_advection_3d(self):
        expr = ['x**2+sin(y)', 'exp(y)+cos(z)', 'tan(z+x)']
        wind = ['tan(y)', 'exp(x)', 'x*exp(z)']
        expect_advection = ['+(tan(y))*(2*x)+(exp(x))*(cos(y))+0',
                            '+0+(exp(x))*(exp(y))+(x*exp(z))*(-sin(z))',
                            '+(tan(y))*(tan(x + z)**2 + 1)'
                            '+0+(x*exp(z))*(tan(x + z)**2 + 1)']
        assert expect_advection == exact_sol.advection(expr, wind)


class ExactSolTests_curl(unittest.TestCase):

    def test_curl_1d(self):
        expr = ['x**2+sin(y)']
        expect_curl = ['cos(y)', '-2*x']
        assert expect_curl == exact_sol.curl(expr)

    def test_curl_2d(self):
        expr = ['x**2+sin(y)', 'exp(y)+cos(z)']
        expect_curl = ['-cos(y)']
        assert expect_curl == exact_sol.curl(expr)

    def test_curl_3d(self):
        expr = ['x**2+sin(y)', 'exp(y)+cos(z)', 'tan(z+x)']
        expect_curl = ['sin(z)', '-tan(x + z)**2 - 1', '-cos(y)']
        assert expect_curl == exact_sol.curl(expr)
