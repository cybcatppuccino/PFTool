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
QP = sympy.expand(x1**5 + x2**5 + x3**5 + x4**5 + x5**5 - 5 * z*x1*x2*x3*x4*x5)
FP = sympy.expand(x1**4 + x2**4 + x3**4 + x4**4 - 4 * z*x1*x2*x3*x4)

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
        outlst = [0 for num in range(len(self.varlist))]
        if inpoly == 0:
            return (outlst, inpoly)
        gbred = sympy.reduced(inpoly, self.gb, gens=self.varlist)
        outpoly = gbred[1]
        for num in range(len(self.gb)):
            for varnum in range(len(self.varlist)):
                outlst[varnum] += gbred[0][num] * self.gbinjac[num][varnum]
        return (outlst, outpoly)
        
    def list_reduced(self, inlst):
        outlst = inlst.copy()
        for deg in range(len(outlst), 0, -1):
            outlst[deg - 1] = sympy.simplify(outlst[deg - 1])
            rs = self.reduced(outlst[deg - 1])
            outlst[deg - 1] = sympy.simplify(rs[1])
            if deg > 1:
                outlst[deg - 2] += sum(sympy.diff(rs[0][num], self.varlist[num]) for num in range(len(self.varlist))) / sympy.Integer(deg - 1)
        while outlst[-1] == 0 and (len(outlst) > 0):
            outlst = outlst[:-1]
        return outlst
    
    def t_derivative(self, inlst):
        outlst = [0 for num in range(len(inlst) + 1)]
        for num in range(len(inlst)):
            outlst[num] += z * sympy.diff(inlst[num], z)
            outlst[num + 1] -= z * (num + 1) * inlst[num] * sympy.diff(self.eqn, z)
        while outlst[-1] == 0 and (len(outlst) > 0):
            outlst = outlst[:-1]
        return self.list_reduced(outlst)
    
    def d_derivative(self, inlst):
        outlst = [0 for num in range(len(inlst) + 1)]
        for num in range(len(inlst)):
            outlst[num] += sympy.diff(inlst[num], z)
            outlst[num + 1] -= (num + 1) * inlst[num] * sympy.diff(self.eqn, z)
        while outlst[-1] == 0 and (len(outlst) > 0):
            outlst = outlst[:-1]
        return self.list_reduced(outlst)

d = PFTool.d
def linear_sol_d(f, g):
    gg = HP(g)
    dlist = [[f]]
    while len(dlist[-1]) == len(dlist):
        dlist.append(gg.d_derivative(dlist[-1]))
    outeqn = d ** (len(dlist) - 1)
    for num in range(len(dlist) - 2, -1, -1):
        if type(dlist[num][num]) == type(1):
            dlist[num][num] = sympy.Integer(dlist[num][num])
        outeqn -= (dlist[-1][num] / dlist[num][num]) * (d ** num)
        for num2 in range(num):
            dlist[-1][num2] -= (dlist[-1][num] / dlist[num][num]) * dlist[num][num2]
    return sympy.simplify(outeqn)

t = PFTool.t
def linear_sol_t(f, g):
    gg = HP(g)
    tlist = [[f]]
    while len(tlist[-1]) == len(tlist):
        tlist.append(gg.t_derivative(tlist[-1]))
    outeqn = t ** (len(tlist) - 1)
    for num in range(len(tlist) - 2, -1, -1):
        if type(tlist[num][num]) == type(1):
            tlist[num][num] = sympy.Integer(tlist[num][num])
        outeqn -= (tlist[-1][num] / tlist[num][num]) * (t ** num)
        for num2 in range(num):
            tlist[-1][num2] -= (tlist[-1][num] / tlist[num][num]) * tlist[num][num2]
    return sympy.simplify(outeqn)

print(linear_sol_d(1, LP))
print(linear_sol_t(1, LP))