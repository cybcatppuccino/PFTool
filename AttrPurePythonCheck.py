'''
(2, '', '4.31', 'This is operator "4.31" from ...') (64*z - 1)**2*(108*z - 1)**2
(3, '60', '4.26', 'Sporadic YY-Operator') (16*z - 1)**2*(108*z - 1)**2
'''

import AESZ
import PFTool
import sympy
import mpmath
from LLL import reduction

#num = 214
#num = 89
#num = 359
#num = 439
num = 190
#num = 175
#num = 3
#num = 371
pfo = AESZ.AESZ(num).pfo
#z0 = sympy.Integer(-1)/20
#z0 = -sympy.Integer(1)/400
#z0 = sympy.Integer(1)/36
#z0 = sympy.Integer(1)/512
z0 = -sympy.Integer(1)
#z0 = -sympy.Integer(1)/7
#z0 = sympy.Integer(1)/8
#z0 = -sympy.Integer(1)/768


print(AESZ.AESZ(num), pfo.calclocalexp())
m, _, r1, r2, r3 = pfo.eval_attr_LLL(z0, pr=True)
print(num, z0, PFTool.mptomma(m), r1, r2, r3)
