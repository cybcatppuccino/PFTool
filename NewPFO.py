import sympy
from PFTool import PFO as PFO
from GetPFO import stupid_linear_sol_t as sol

x1 = sympy.Symbol('x1')
x2 = sympy.Symbol('x2')
x3 = sympy.Symbol('x3')
x4 = sympy.Symbol('x4')
z = sympy.Symbol('z')
t = sympy.Symbol('t')

#C0 = sympy.expand((x1+x2+x3)*(x1+x2)*(x1+x2-x3)*z+x1*x2*x3) # 1/4, 3/4
#C1 = sympy.expand((x1**3+x2*(x2+x3)*(x2-x3))*z+x1*x2*x3) # 1/6, 5/6
#C = sympy.expand((1-4*z)*x1*x2*x3-z*(x1+x2)*(x1*x2+x3*x3))
#C = sympy.expand((z-1728)*(x2*x2*x3-4*x1*x1*x1)+27*z*x1*x3*x3+27*z*x3*x3*x3)
C = sympy.expand((x1**3+x2**3+x3**3)*z-3*x1*x2*x3)
r = sol(1, C)
print("r=",r)
p = PFO(r.as_numer_denom()[0])
print("p=",p)
print(p.all_sol(30))
print(p.calclocalexp())
'''
C = sympy.expand((x1+x2+x3)*(x1*x2+x2*x3+x3*x1)-z*x1*x2*x3)
r = sol(1, C)
print("r=",r)
q = PFO(r.as_numer_denom()[0])
#print("p=",q)

print(q.all_sol(30))
print(q.calclocalexp())
'''