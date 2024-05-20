import sympy
import fractions

z = sympy.Symbol('z')
print(1)
f1 = sympy.fps((135*z**4 + 292*z**3 + 178*z**2 + 20*z - 1)/(15*z + 43))
print(f1.truncate(8).removeO())

