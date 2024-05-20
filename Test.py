from PFTool import *
from AESZ import *
from sympy import *

x = Symbol('x')
y = Symbol('y')
z = Symbol('z')
gb = groebner([x ** 2 + y, y ** 2 + 3 * x + z], [x,y], domain='QQ(z)')
gb.a
print(gb)
print(list(gb))

