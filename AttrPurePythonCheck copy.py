'''
(2, '', '4.31', 'This is operator "4.31" from ...') (64*z - 1)**2*(108*z - 1)**2
(3, '60', '4.26', 'Sporadic YY-Operator') (16*z - 1)**2*(108*z - 1)**2
'''

import AESZ
import PFTool
import sympy
import mpmath
from LLL import reduction

num = 414
#num = 175
#num = 3
#num = 371
pfo = AESZ.AESZ(num).pfo
print(pfo, pfo.calclocalexp())
m = pfo.monodromy_only_pt(1/64, pr=True)
print(PFTool.mptomma(m))
