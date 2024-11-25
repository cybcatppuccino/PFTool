'''
(2, '', '4.31', 'This is operator "4.31" from ...') (64*z - 1)**2*(108*z - 1)**2
(3, '60', '4.26', 'Sporadic YY-Operator') (16*z - 1)**2*(108*z - 1)**2
'''

import AESZ
import PFTool
import sympy
import mpmath
from LLL import reduction

num = 227
pfo = AESZ.AESZ(num).pfo
z0 = sympy.Integer(1)/2

print(pfo.calclocalexp())

def LLLtest(mat):
    ls = list((mat ** -1) * mpmath.matrix([[1], [0], [0], [0]]))
    tpj = 2 * mpmath.pi * mpmath.j
    vals = [mpmath.zeta(3), tpj ** 3, tpj ** 2, tpj ** 1]
    tol = 10**(12 - mpmath.mp.dps)
    ls2 = [mpmath.re(ls[0]*vals[0]), mpmath.re(ls[0]*vals[1]), mpmath.re(ls[1]*vals[2]), mpmath.re(ls[2]*vals[3]), mpmath.re(ls[3]), \
        mpmath.im(ls[0]*vals[0]), mpmath.im(ls[0]*vals[1]), mpmath.im(ls[1]*vals[2]), mpmath.im(ls[2]*vals[3]), mpmath.im(ls[3])]
    ls3 = []
    for _ in ls2:
        if abs(_) > tol:
            ls3.append(_)
    l = len(ls3)
    llllst = [[int(a==b) + int(ls3[a] * 10**(mpmath.mp.dps - 12)) * int(b==l) for b in range(l+1)] for a in range(l)]
    return ls2, reduction(llllst)

def LLLtest2(mat):
    ls = list((mat ** -1) * mpmath.matrix([[1], [0], [0], [0]]))
    tpj = 2 * mpmath.pi * mpmath.j
    vals = [mpmath.zeta(3), tpj ** 3, tpj ** 2, tpj ** 1]
    ls2 = [mpmath.re(ls[0]*vals[0]), mpmath.re(ls[0]*vals[1]), mpmath.re(ls[1]*vals[2]), mpmath.re(ls[2]*vals[3]), mpmath.re(ls[3]), \
        mpmath.im(ls[0]*vals[0]), mpmath.im(ls[0]*vals[1]), mpmath.im(ls[1]*vals[2]), mpmath.im(ls[2]*vals[3]), mpmath.im(ls[3])]
    l = len(ls2) // 2
    llllst = [[int(a==b) + int(ls2[a] * 10**(mpmath.mp.dps - 12)) * int(b==l) + int(ls2[a + l] * 10**(mpmath.mp.dps - 12)) * int(b==l+1) for b in range(l+2)] for a in range(l)]
    return ls2, reduction(llllst)

mat = pfo.eval_Wronskian0(z0, pr=True)
rst1, rst2 = LLLtest(mat)
print(num, z0, rst2)
rst1, rst2 = LLLtest2(mat)
print(num, z0, rst2)