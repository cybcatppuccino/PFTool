
import AESZ
import PFTool
import sympy
import mpmath
from LLL import reduction

num = 371
pfo = AESZ.AESZ(num).pfo
coefflist = list(map(int, sympy.Poly(pfo.discriminant).all_coeffs()))
print(coefflist)
mroot = min(mpmath.polyroots(coefflist, maxsteps=10000), key=abs)
m = pfo.monodromy_only_pt(mroot, pr=True)
bc = PFTool.BASECHANGE
const = mpmath.zeta(3) / ((2 * mpmath.pi)**3)
m1 = bc * m * (bc ** -1) - mpmath.eye(4)
print(-mpmath.re(m1[3, 0]), -mpmath.re(m1[1, 0])*24, mpmath.im(m1[0, 0]) / const)
