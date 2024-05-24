import math
import sympy
import fractions
from numba import jit, njit
import numba

def foo(a, b):
    for _ in range(1000000):
        fractions.Fraction(a, b) + 1
    return fractions.Fraction(a, b) + 1


def foo2(a, b):
    for _ in range(1000000):
        fractions.Fraction(a, b) + 1
    return fractions.Fraction(a, b) + 1


def int_to_p_adic(inint, inprime):
    if not inint:
        return (0, math.inf)
    outtuple = (inint, 0)
    while True:
        rs = divmod(outtuple[0], inprime)
        if rs[1]:
            return outtuple
        else:
            outtuple = (rs[0], outtuple[1] + 1)
            
def eGCD(a, b):
    prevx, x = 1, 0
    prevy, y = 0, 1
    while b:
        rs = divmod(a, b)
        x, prevx = prevx - rs[0]*x, x
        y, prevy = prevy - rs[0]*y, y
        a, b = b, rs[1]
    if a >= 0:
        return a, (prevx, prevy)
    else:
        return -a, (-prevx, -prevy)

def rat_to_p_adic(inint1, inint2, inprime, inacc):
    rs1 = int_to_p_adic(inint1, inprime)
    rs2 = int_to_p_adic(inint2, inprime)
    po = inprime ** inacc
    return ((rs1[0] * eGCD(rs2[0], po)[1][0]) % po, rs1[1] - rs2[1])

def to_p_adic(term, inprime, inacc):
    if type(term) == type(1):
        return int_to_p_adic(term, inprime)
    elif type(term) == type(fractions.Fraction(1, 2)):
        return rat_to_p_adic(int(term.numerator), int(term.denominator), inprime, inacc)
    elif term.is_Integer:
        return int_to_p_adic(int(term), inprime)
    elif term.is_Rational:
        return rat_to_p_adic(int(term.numerator()), int(term.denominator()), inprime, inacc)

# The p-adic number Class
class PN:
    def __init__(self, inprime, inacc):
        self.p = inprime
        self.acc = inacc
        self.ppow = inprime ** inacc
    def setint(self, inint):
        rs = int_to_p_adic(inint, self.p)
        self.deg = rs[1]
        self.num = rs[0] % self.ppow
        return self
    def setfrac(self, innum, inden):
        rs1 = int_to_p_adic(innum, self.p)
        rs2 = int_to_p_adic(inden, self.p)
        self.deg = rs1[1] - rs2[1]
        self.num = (rs1[0] * eGCD(rs2[0], self.ppow)[1][0]) % self.ppow
        return self
    def setval(self, term):
        if type(term) == type(1):
            return self.setint(term)
        elif type(term) == type(fractions.Fraction(1, 2)):
            return self.setfrac(int(term.numerator), int(term.denominator))
        elif term.is_Integer:
            return self.setint(int(term))
        elif term.is_Rational:
            return self.setfrac(int(term.numerator()), int(term.denominator()))
        else:
            raise Exception("Type!")
    
    def __int__(self):
        if self.deg >= self.acc:
            return 0
        elif self.deg >= 0:
            return (self.num * (self.p ** self.deg)) % self.ppow
        else:
            return 0
    
    def __str__(self):
        if self.deg == math.inf:
            return '0'
        if self.deg == 0:
            aux = ''
        elif self.deg == 1:
            aux = 'p'
        else:
            aux = 'p^' + str(self.deg)
        return str(self.num) + aux
    
    def __repr__(self):
        b = 'PN(' + str(self.p) + ', ' + str(self.acc)
        if self.deg == math.inf:
            return b + ').setint(0)'
        elif self.deg == 0:
            return b + ').setint(' + str(self.num + (self.acc == 0)) + ')'
        elif self.deg == 1:
            return b + ').setint(' + str(self.num + (self.acc == 0)) + ' * ' + str(self.p) + ')'
        elif self.deg > 0:
            return b + ').setint(' + str(self.num + (self.acc == 0)) + ' * ' + str(self.p) + '**' + str(self.deg) + ')'
        else:
            return b + ').setfrac(' + str(self.num + (self.acc == 0)) + ', ' + str(self.p) + '**' + str(-self.deg) + ')'
    
    
    def __mul__(self, term):
        if type(term) == type(self):
            outPN = PN(self.p, min(self.acc, term.acc))
            outPN.deg = self.deg + term.deg
            outPN.num = (self.num * term.num) % outPN.ppow
        else:
            rs = to_p_adic(term, self.p, self.acc)
            outPN = PN(self.p, self.acc)
            outPN.deg = self.deg + rs[1]
            outPN.num = (self.num * rs[0]) % outPN.ppow
        return outPN
            
    def __rmul__(self, term):
        return self * term
    
    
    def __truediv__(self, term):
        if type(term) == type(self):
            outPN = PN(self.p, min(self.acc, term.acc))
            outPN.deg = self.deg - term.deg
            outPN.num = (self.num * eGCD(term.num, self.ppow)[1][0]) % outPN.ppow
        else:
            rs = to_p_adic(term, self.p, self.acc)
            outPN = PN(self.p, self.acc)
            outPN.deg = self.deg - rs[1]
            outPN.num = (self.num * eGCD(rs[0], self.ppow)[1][0]) % outPN.ppow
        return outPN
    
    
    def __rtruediv__(self, term):
        if type(term) == type(self):
            outPN = PN(self.p, self.acc)
            outPN.deg = term.deg - self.deg
            outPN.num = (term.num * eGCD(self.num, self.ppow)[1][0]) % self.ppow
        else:
            rs = to_p_adic(term, self.p, self.acc)
            outPN = PN(self.p, self.acc)
            outPN.deg = - self.deg + rs[1]
            outPN.num = (rs[0] * eGCD(self.num, self.ppow)[1][0]) % self.ppow
        return outPN
    
    def copy(self):
        outPN = PN(self.p, self.acc)
        outPN.num = self.num
        outPN.deg = self.deg
        return outPN
    
    def result_de_acr(self):
        rs = int_to_p_adic(self.num, self.p)
        self.deg += rs[1]
        self.acc -= rs[1]
        self.ppow = self.p ** self.acc
        self.num = rs[0] % self.ppow
    
    
    def self_de_acr(self, absacr):
        self.acc = absacr - self.deg
        self.ppow = self.p ** self.acc
        self.num = self.num % self.ppow
    
    def __add__(self, term):
        if type(term) == type(self):
            if term.p != self.p:
                raise Exception("Different Prime!")
            if self.deg <= term.deg:
                absac1 = self.deg + self.acc
                if term.deg >= absac1:
                    return self.copy()
                absac2 = term.deg + term.acc
                outPN = self.copy()
                if absac2 < absac1:
                    outPN.self_de_acr(absac2)
                outPN.num += term.num * self.p ** (term.deg - self.deg)
                outPN.result_de_acr()
            else:
                absac1 = term.deg + term.acc
                if self.deg >= absac1:
                    return term.copy()
                absac2 = self.deg + self.acc
                outPN = term.copy()
                if absac2 < absac1:
                    outPN.self_de_acr(absac2)
                outPN.num += self.num * term.p ** (self.deg - term.deg)
                outPN.result_de_acr()
        else:
            rs = to_p_adic(term, self.p, self.acc)
            if self.deg <= rs[1]:
                absac1 = self.deg + self.acc
                if rs[1] >= absac1:
                    return self.copy()
                outPN = self.copy()
                outPN.num += rs[0] * self.p ** (rs[1] - self.deg)
                outPN.result_de_acr()
            else:
                outPN = PN(self.p, self.acc)
                if self.deg >= rs[1] + self.acc:
                    outPN.num = rs[0] % self.ppow
                else:
                    outPN.num = (rs[0] + self.num * self.p ** (self.deg - rs[1])) % self.ppow
                outPN.deg = rs[1]
        return outPN
    
    def __radd__(self, term):
        return self + term
    
    def __neg__(self):
        outPN = self.copy()
        outPN.num = (-outPN.num) % outPN.ppow
        return outPN
    
    def __sub__(self, term):
        return self + (-term)
    
    def __rsub__(self, term):
        return (-self) + term
    
    def __bool__(self):
        return self.deg != math.inf
    
    
    def __pow__(self, num):
        num = int(num)
        if num == 0:
            return PN(self.p, self.acc).setint(1)
        elif num == 1:
            return self.copy()
        elif num > 1:
            if num % 2 == 0:
                a = self ** (num / 2)
                return a * a
            else:
                a = self ** ((num-1) / 2)
                return a * a * self
        elif num < 0:
            a = 1 / self
            return a ** (-num)

def Teichmuller(a, p, n):
    a = a % p
    if a == 1:
        return PN(p, n).setint(1)
    elif a == p-1:
        return PN(p, n).setint(-1)
    else:
        outPN = PN(p, n).setint(a)
        for _ in range(n - 1):
            outPN **= p
        return outPN

if __name__ == "__main__":
    a = PN(3,5).setint(1)
    # a = 1
    print("started")
    for _ in range(100000):
        a += a
    print("finished")
