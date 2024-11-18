import AESZ
import sympy

'''
Attractor Point
level N = 27
infty index a1, a2, a3, a4 = 1/3, 1/3, 2/3, 2/3
conifold pt = 1/3^6 = 1/729
mirror M = X_{3,3}(1^6)

kappa = 9
c2.D = 54
Euler chi(M) = -144

AESZ: (207, '4', '1.4', 'A-incarnation: X(3,3) in P^5.')
'''

X = AESZ.AESZ(207)
op = X.pfo
print(op)
# print(op.calclocalexp())
# (1/729, {2: 1, 1: 2, 0: 1}), (0, {0: 4}), (zoo, {2/3: 2, 1/3: 2})

print(op.all_sol(6))
# Holomorphic Solution: (3n)!^2 / n!^6
at = -sympy.Integer(1) / 5832
opat = op.translation(at)
print(opat)