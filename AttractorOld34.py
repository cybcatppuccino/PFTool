import AESZ
import PFTool
import sympy
import sys
sys.set_int_max_str_digits(100000)

'''
Attractor Point

z* = -1/768

level N = ?
conifold pt = 1/2^8 = 1/256
mirror M = ?

kappa = ?
c2.D = ?
Euler chi(M) = ?

AESZ: (371, '111', '2.17', 'B-Incarnations:')
'''

X = AESZ.AESZ(379)
op = X.pfo
print(op)
print(op.calclocalexp())

op1 = op.translation_inf()
op2 = op1.translation(864)
op3 = op2.translation_inf()
newtform = sympy.simplify(op3.tform * (864*PFTool.z+1))
op4 = PFTool.PFO(newtform)
print(op4)


solinfo = str(op4.all_sol(800)[0])
solinfo = solinfo.replace('[', '{')
solinfo = solinfo.replace(']', '}')

with open('solinnew.txt', 'w') as f:
    f.write('sol=' + solinfo + ';')

