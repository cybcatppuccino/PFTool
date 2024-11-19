import AESZ
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

X = AESZ.AESZ(371)
op = X.pfo
print(op)
print(op.calclocalexp())
# (1/256, {1: 1, 1/2: 2, 0: 1}), (0, {0: 4}), (zoo, {3/2: 2, 1/2: 2})]


solinfo = op.all_sol(800)[0]

with open('solinnew.txt', 'w') as f:
    f.write(str(solinfo))
