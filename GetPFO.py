import sympy
import PFTool

# sympy.Symbols
x1 = sympy.Symbol('x1')
x2 = sympy.Symbol('x2')
x3 = sympy.Symbol('x3')
x4 = sympy.Symbol('x4')
x5 = sympy.Symbol('x5')

# extra variables to represent the Groebner basis into the original one
y1 = sympy.Symbol('y1')
y2 = sympy.Symbol('y2')
y3 = sympy.Symbol('y3')
y4 = sympy.Symbol('y4')
y5 = sympy.Symbol('y5')

# the family variable
z = sympy.Symbol('z')

# Legendre Homogeneous Polynomial
LP = sympy.expand(x3 * x2**2 + x1*(x1 - x3)*(x1 - z * x3))
QP = sympy.expand(x1**5 + x2**5 + x3**5 + x4**5 + x5**5 + z*x1*x2*x3*x4*x5)
FP = sympy.expand(x1**4 + x2**4 + x3**4 + x4**4 + z*x1*x2*x3*x4)

# We always consider homogeneous polynomials

'''
def homo_poly_coeff_list(poly, varlist, deg):
    if len(varlist) == 1:
        return [poly.coeff(varlist[0], deg)]
    else:
        var = varlist[0]
        varlist2 = varlist[1:]
        outlst = []
        for num in range(deg, -1, -1):
            outlst += homo_poly_coeff_list(poly.coeff(var, num), varlist2, deg - num)
        return outlst
print(homo_poly_coeff_list(LP, [x1,x2,x3], 3))

def jac_ideal(poly, varlist):
    return [sympy.simplify(poly.diff(var)) for var in varlist]
'''

class HP:
    def __init__(self, ineqn):
        self.eqn = sympy.expand(ineqn)
        self.varlist = []
        self.auxvarlist = []
        for num in range(5):
            if "x" + str(num + 1) in str(self.eqn):
                self.varlist.append([x1, x2, x3, x4, x5][num])
                self.auxvarlist.append([y1, y2, y3, y4, y5][num])
        self.jac = [sympy.simplify(self.eqn.diff(self.varlist[num])) for num in range(len(self.varlist))]
        self.auxjac = [self.jac[num] - self.auxvarlist[num] for num in range(len(self.varlist))]
        self.auxgb = sympy.groebner(self.auxjac, self.varlist + self.auxvarlist, domain='QQ(z)', order='grlex')
        
        self.gb = []
        self.gbinjac = []
        jacsubs = [(self.auxvarlist[num], self.jac[num]) for num in range(len(self.varlist))]
        for term0 in self.auxgb:
            term = term0
            outlst = []
            for var in self.auxvarlist:
                const = term.coeff(var, 0)
                rest = term - const
                term = const
                outlst.append(sympy.cancel(rest / var).subs(jacsubs))
            self.gb.append(-term)
            self.gbinjac.append(outlst)
            
    def reduced(self, inpoly):
        
        
        
l = HP(FP)
