import sympy
import PFTool
import MyGroebner

# sympy.Symbols
x1 = sympy.Symbol('x1')
x2 = sympy.Symbol('x2')
x3 = sympy.Symbol('x3')
x4 = sympy.Symbol('x4')
x5 = sympy.Symbol('x5')
x6 = sympy.Symbol('x6')

# the family variable
z = sympy.Symbol('z')

# Legendre Homogeneous Polynomial
LP2 = sympy.expand(x3 * x2**2 + x1*(z * x1 + x3)*(x1 - z * x3))
LP = sympy.expand(x3 * x2**2 + x1*(x1 - x3)*(x1 - z * x3))
QP = sympy.expand(x1**5 + x2**5 + x3**5 + x4**5 + x5**5 - (1 / z)*x1*x2*x3*x4*x5)
CP3 = sympy.expand(((x1+x2+x3)*(1/x1+1/x2+1/x3)*z-1)*x1*x2*x3)
CP4 = sympy.expand(((x1+x2+x3)*(1/x1+1/x2+1/x3)-z)*x1*x2*x3)
QP2 = sympy.expand(((x1+x2+x3+x4+x5)*(1/x1+1/x2+1/x3+1/x4+1/x5)-z)*x1*x2*x3*x4*x5)
FP2 = sympy.expand(x1**4 + x2**4 + x3**4 + x4**4 - (1 / z)*(x3**2*x4**2 + x1*x2*x3*x4))
FP = sympy.expand(x1**4 + x2**4 + x3**4 + x4**4 - (1 / z)*x1*x2*x3*x4)
CP2 = sympy.expand(x1**3 + x2**3 + x3**3 - (1 / z)*(x1**2 * x2 + x1*x2*x3))
CP = sympy.expand(x1**3 + x2**3 + x3**3 - (1 / z)*x1*x2*x3)

C0 = sympy.expand(x1*(x1-x3)*(x2-x3) + z*(x1-x2)*x2*x3)



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
        # self.auxvarlist = []
        for num in range(6):
            if "x" + str(num + 1) in str(self.eqn):
                self.varlist.append([x1, x2, x3, x4, x5, x6][num])
                # self.auxvarlist.append([y1, y2, y3, y4, y5][num])
        
        print("Preparation done!")
        
        self.jac = [sympy.simplify(self.eqn.diff(self.varlist[num])) for num in range(len(self.varlist))]
        # self.auxjac = [self.jac[num] - self.auxvarlist[num] for num in range(len(self.varlist))]
        
        print(str(self.jac).replace("**", "^"), self.varlist)
        
        # print(MyGroebner.GroebnerBasis(self.jac, self.varlist, domain='QQ(z)', order='grlex'))
        # self.auxgb = MyGroebner.GroebnerBasis(self.auxjac, self.varlist + self.auxvarlist, domain='QQ(z)', order='grlex')
        
        self.gb, self.gbinjac = MyGroebner.GroebnerBasis(self.jac, self.varlist, domain='QQ(z)', order='grlex')
        
        print("extendedGB done! Length = ", len(self.gb))
            
    def reduced(self, inpoly):
        print("Reducing ", inpoly)
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
        print("Reducing list", str(inlst)[:50])
        outlst = inlst.copy()
        for deg in range(len(outlst), 0, -1):
            # outlst[deg - 1] = sympy.simplify(outlst[deg - 1])
            rs = self.reduced(outlst[deg - 1])
            # outlst[deg - 1] = sympy.simplify(rs[1])
            outlst[deg - 1] = rs[1]
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
def stupid_linear_sol_d(f, g):
    gg = HP(g)
    dlist = [[f]]
    while len(dlist[-1]) == len(dlist):
        dlist.append(gg.d_derivative(dlist[-1]))
    print(dlist)
    outeqn = d ** (len(dlist) - 1)
    for num in range(len(dlist) - 2, -1, -1):
        if type(dlist[num][num]) == type(1):
            dlist[num][num] = sympy.Integer(dlist[num][num])
        outeqn -= (dlist[-1][num] / dlist[num][num]) * (d ** num)
        for num2 in range(num):
            dlist[-1][num2] -= (dlist[-1][num] / dlist[num][num]) * dlist[num][num2]
    return sympy.simplify(outeqn)

t = PFTool.t
def gen_t_list(f, g, n):
    gg = HP(g)
    tlist = [[f]]
    for num in range(n):
        tlist.append(gg.t_derivative(tlist[-1]))
    print(tlist)
    return tlist

def stupid_linear_sol_t(f, g):
    gg = HP(g)
    tlist = [[f]]
    while len(tlist[-1]) == len(tlist):
        tlist.append(gg.t_derivative(tlist[-1]))
    print(tlist)
    outeqn = t ** (len(tlist) - 1)
    for num in range(len(tlist) - 2, -1, -1):
        if type(tlist[num][num]) == type(1):
            tlist[num][num] = sympy.Integer(tlist[num][num])
        outeqn -= (tlist[-1][num] / tlist[num][num]) * (t ** num)
        for num2 in range(num):
            tlist[-1][num2] -= (tlist[-1][num] / tlist[num][num]) * tlist[num][num2]
    return sympy.simplify(outeqn)


# gen_t_list(x1, x3 * LP, 2)
# gb = MyGroebner.GroebnerBasis([x1**3, x2**3], [x1, x2], domain='QQ(z)', order='grlex')

C3 = sympy.expand((x1**3+x2**3+x3**3)*z-3*x1*x2*x3)
C3b = sympy.expand((x1*x1*x2+x2*x2*x3+x3*x3*x1)*z-x1*x2*x3)
C6 = sympy.expand((x1+x2+x3)*(x1*x2+x2*x3+x3*x1) - z*x1*x2*x3)
C2 = sympy.expand((x1+x2)*(x1*x2+x3*x3)-16*z*x1*x2*x3)
C5 = sympy.expand((x1+x2)*(x1+x3)*(x1+x2+x3)+z*x1*x2*x3)
C6b = sympy.expand((x1+x2+x3)*(x1*x2+x2*x3+x3*x3) + 8*z*x1*x2*x3)
C9 = sympy.expand(x1**3+x2*x3*(x2+x3)+(27*z-3)*x1*x2*x3)

Ct = sympy.expand((x1*x1+1)*(x2*x2+1)*(x3*x3+1)-z)
#rest = stupid_linear_sol_d(1, Ct)

#res2 = stupid_linear_sol_d(1, C2)
#res3 = stupid_linear_sol_d(1, C3)
#res3b = stupid_linear_sol_d(1, C3b)
#res6 = stupid_linear_sol_d(1, C6)
#res5 = stupid_linear_sol_d(1, C5)
#res6b = stupid_linear_sol_d(1, C6b)
#res9 = stupid_linear_sol_d(1, C9)

p3 = PFTool.PFO("(d**2*z*(z**3 - 1) + d*(4*z**3 - 1) + 2*z**2)")
p3b = PFTool.PFO("(27*d**2*z**4 - d**2*z + 108*d*z**3 - d + 54*z**2)")
p6 = PFTool.PFO("(d**2*z**3 - 10*d**2*z**2 + 9*d**2*z + 3*d*z**2 - 20*d*z + 9*d + z - 3)")
p2 = PFTool.PFO("(16*d**2*z**3 - d**2*z + 48*d*z**2 - d + 16*z)")
p5 = PFTool.PFO("(d**2*z**3 + 11*d**2*z**2 - d**2*z + 3*d*z**2 + 22*d*z - d + z + 3)")
p6b = PFTool.PFO("(8*d**2*z**3 + 7*d**2*z**2 - d**2*z + 24*d*z**2 + 14*d*z - d + 8*z + 2)")
p9 = PFTool.PFO("(27*d**2*z**3 - 9*d**2*z**2 + d**2*z + 81*d*z**2 - 18*d*z + d + 27*z - 3)")

ptest = PFTool.PFO("(-z^2+34*z^3-z^4)*d**3 + (-3*z+153*z^2-6*z^3)*d**2 + (-1+112*z-7*z^2)*d + (5-z)")
#print(p3.all_sol(20)[0])
#pa = PFTool.PFO("(d**2*z**3 - 10*d**2*z**2 + 9*d**2*z + 3*d*z**2 - 20*d*z + 9*d + z - 3)")
#pb = PFTool.PFO("(d**2*z**5 - 28*d**2*z**4 + 134*d**2*z**3 + 100*d**2*z**2 - 63*d**2*z + 3*d*z**4 - 54*d*z**3 + 76*d*z**2 + 326*d*z - 63*d + z**3 - 7*z**2 - 33*z + 39)")

# print(stupid_linear_sol_d(1, QP))
# print(stupid_linear_sol_t(1, FP))
# gen_t_list(1, FP2, 4)
# print(len(sympy.GroebnerBasis([x2*x3*(2*x1*x2*x3*z - 2*x1*x2*x4*z - 2*x1*x3*x4*z + 2*x1*x4**2*z - x2*x3*x4*z + x2*x4**2*z + x3*x4**2*z - x4**3*z + x4**3), x1*x3*(2*x1*x2*x3*z - 2*x1*x2*x4*z - x1*x3*x4*z + x1*x4**2*z - 2*x2*x3*x4*z + 2*x2*x4**2*z + x3*x4**2*z - x4**3*z + x4**3), 2*x1**2*x2**2*x3*z - x1**2*x2**2*x4*z - 2*x1**2*x2*x3*x4*z + x1**2*x2*x4**2*z - 2*x1*x2**2*x3*x4*z + x1*x2**2*x4**2*z + 2*x1*x2*x3*x4**2*z - x1*x2*x4**3*z + x1*x2*x4**3 - x4**5, -x1**2*x2**2*x3*z - x1**2*x2*x3**2*z + 2*x1**2*x2*x3*x4*z - x1*x2**2*x3**2*z + 2*x1*x2**2*x3*x4*z + 2*x1*x2*x3**2*x4*z - 3*x1*x2*x3*x4**2*z + 3*x1*x2*x3*x4**2 - 5*x3*x4**4 + 6*x4**5], [x1, x2, x3, x4], domain='QQ(z)', order='grlex')))
# print(stupid_linear_sol_d(x4**2, x4**6 - x4**3 * (x4**2 - x1*x2) * x3 - z * x1*x2*x3 * (x4-x1)*(x4-x2)*(x4-x3)))
