'''
(2, '', '4.31', 'This is operator "4.31" from ...') (64*z - 1)**2*(108*z - 1)**2
(3, '60', '4.26', 'Sporadic YY-Operator') (16*z - 1)**2*(108*z - 1)**2
'''

import AESZ
import PFTool
import sympy
import mpmath
from LLL import reduction

# PFO100STR = "t^4 - z * (73*t^4+98*t^3+77*t^2+28*t+4) + z^2 * (520*t^4-1040*t^3-2904*t^2-2048*t-480) + 2^6 * z^3 * (65*t^4+390*t^3+417*t^2+180*t+28) - 2^9 * z^4 * (73*t^4+194*t^3+221*t^2+124*t+28)+2^15 * z^5 * (t+1)^4"

# num = 214
# num = 89
# num = 359
# num = 439
# num = 190
num = 175
# num = 3
# num = 371
# num = 175

pfo = AESZ.AESZ(num).pfo

# z0 = sympy.Integer(-1)/20
# z0 = -sympy.Integer(1)/400
# z0 = sympy.Integer(1)/36
# z0 = sympy.Integer(1)/512
# z0 = -sympy.Integer(1)
z0 = sympy.Integer(1)/3
# z0 = sympy.Integer(1)/8
# z0 = -sympy.Integer(1)/768
# z0 = sympy.Integer(33) - 8 * sympy.sqrt(17)

val = pfo.attrval_MUM(200, z0)
idt = mpmath.identify(val, [], tol=10**-30)
print(val, idt)

# pfo = PFTool.PFO(PFO100STR)
# z0 = sympy.Integer(1)/8

# print(AESZ.AESZ(num), pfo.calclocalexp())
# print(pfo.calclocalexp())
m, _, r1, r2, r3 = pfo.eval_attr_LLL(z0, pr=True)
print(z0, PFTool.mptomma(m), r1, r2, r3)


