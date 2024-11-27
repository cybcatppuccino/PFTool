'''
(2, '', '4.31', 'This is operator "4.31" from ...') (64*z - 1)**2*(108*z - 1)**2
(3, '60', '4.26', 'Sporadic YY-Operator') (16*z - 1)**2*(108*z - 1)**2
'''

import AESZ
import PFTool
import sympy
import mpmath
from LLL import reduction

num = 3
pfo = AESZ.AESZ(num).pfo
z0 = sympy.Integer(1)/8 # + 0 * sympy.sqrt(sympy.Integer(17))

print(pfo.calclocalexp())
print(num, z0, pfo.eval_attr_LLL(z0, pr=True))
print(num, z0, pfo.eval_Kp(z0, pr=True))