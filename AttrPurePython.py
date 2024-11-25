import AESZ
import sympy
import mpmath

pfo = AESZ.AESZ(447).pfo
z0 = sympy.Integer(1)/8
#print(pfo.eval_attr_Wronskian_MMA(z0, pr=True))


print(pfo.eval_attr_LLL(z0, pr=True))