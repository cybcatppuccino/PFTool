'''
Picard-Fuchs Operator Tool by cybcat

We only deal with PF Operator in polynomial of z, d = d/dz and t = z * d/dz.
'''

import sympy

# The Symbols
# z: the variable of the ODE
# lgz: log(z), only used in log solutions
# d: the differential operator d/dz
# t: the differential operator z * d/dz

# zdelta: used when considering z to z+zdelta

z = sympy.Symbol('z')
lgz = sympy.Symbol('lgz')
d = sympy.Symbol('d')
t = sympy.Symbol('t')

zdelta = sympy.Symbol('zdelta')

DEG_BOUND = 15
Z_DEG_BOUND = 100
TEST_PFO = "(t ** 4) - z * (t + 1/5) * (t + 2/5) * (t + 3/5) * (t + 4/5)"
LEG = "z * (z-1) * d^2 + (2*z - 1) * d + 1/4"


# We introduce
X = sympy.Symbol('t')
# Because on AESZ List they have (e.g. AESZ300 and AESZ22): 
TEST_LIST300 = [X**4, 240+2560*X+9456*X**2+13792*X**3+5936*X**4, 
1628160+13957120*X+33556480*X**2+21186560*X**3+2293760*X**4,
-3422617600*X**4-12288000000*X**3-9065267200*X**2-2457600000*X-221184000,
1048576000*(5*X+1)*(5*X+2)*(5*X+3)*(5*X+4)]
TEST_LIST22 = [49*X**4, -1085*X**4-2002*X**3-1638*X**2-637*X-98,
-16105*X**4-68044*X**3-102261*X**2-66094*X-15736, 21000*X**4+68712*X**3+72568*X**2+30072*X+3808,
-7440*X**4-20256*X**3-23024*X**2-12896*X-2944, 512*(X+1)**4]


# Generate a list of power of t in d-poly form
def t_power_to_d_list(n):
    lst = [1]
    for _ in range(n):
        lst.append(sympy.expand((lst[-1] * d + sympy.diff(lst[-1], z)) * z))
    return lst

# Generate a list of power of d in t-poly form
def d_power_to_t_list(n):
    lst = [1]
    for _ in range(n):
        lst.append(sympy.expand(lst[-1] * t / z + sympy.diff(lst[-1], z)))
    return lst

# Generate a list of power of d when changing variable from z to 1/z
def d_power_inv_list(n):
    lst = [1]
    for _ in range(n):
        lst.append(sympy.expand((lst[-1] * d + sympy.diff(lst[-1], z)) * (- z**2)))
    return lst

def to_d_form(ineqn):
    ineqn = sympy.expand(ineqn)
    t_deg = DEG_BOUND
    while (ineqn.coeff(t, t_deg) == 0):
        t_deg -= 1
    lst = t_power_to_d_list(t_deg)
    outeqn = 0
    for deg in range(t_deg + 1):
        outeqn += ineqn.coeff(t, deg) * lst[deg]
    return sympy.expand(outeqn)

def to_t_form(ineqn):
    ineqn = sympy.expand(ineqn)
    d_deg = DEG_BOUND
    while (ineqn.coeff(d, d_deg) == 0):
        d_deg -= 1
    lst = d_power_to_t_list(d_deg)
    outeqn = 0
    for deg in range(d_deg + 1):
        outeqn += ineqn.coeff(d, deg) * lst[deg]
    return sympy.expand(outeqn)

# d form inverse when changing variable from z to 1/z
def to_inv(ineqn, indeg):
    ineqn = sympy.expand(ineqn)
    lst = d_power_inv_list(indeg)
    outeqn = 0
    for deg in range(indeg + 1):
        outeqn += ineqn.coeff(d, deg).subs([(z, 1/z)]) * lst[deg]
    return sympy.expand(outeqn)

def to_primitive(ineqn, var, deg):
    ineqn = sympy.expand(ineqn * (z ** DEG_BOUND))
    gcd = ineqn.coeff(var, 0)
    for num in range(1, deg + 1):
        gcd = sympy.gcd(gcd, ineqn.coeff(var, num))
    return sympy.expand(sympy.cancel(ineqn / gcd))

# The Picard Fuchs Operator Class
class PFO:
    def __init__(self, ineqn, pr=False):
        self.pr = pr
        if pr:
            print("Constructing Operator...")
        
        # Considering String Input
        eqn = ineqn
        if type(ineqn) == type("string"):
            eqn = sympy.sympify(ineqn, evaluate=False, rational=True, convert_xor=True)
        elif type(ineqn) == type(["List", 114514]):
            eqn = 0
            for num in range(len(ineqn)):
                eqn += ineqn[num] * (z ** num)
        ineqn = sympy.expand(eqn)
        
        # Standard Forms
        self.dform = to_d_form(ineqn)
        self.tform = to_t_form(ineqn)
        
        # Degree
        self.deg = DEG_BOUND
        while (self.dform.coeff(d, self.deg) == 0):
            self.deg -= 1
        
        # Primitive Forms
        self.primdform = to_primitive(self.dform, d, self.deg)
        self.primtform = to_primitive(self.tform, t, self.deg)
        
        # Indicial Polynomial
        poly = self.primtform.subs([(z, 0)])
        self.localind = sympy.factor(sympy.cancel(poly / poly.coeff(t, self.deg)))
        
        # Discriminant
        self.discriminant = sympy.factor(self.primtform.coeff(t, self.deg))
        
        # Local Exponents (Not calculated at first)
        self.localexp = None
        # Primitive t-form List
        self.primtlist = None
        
        if pr:
            print("Operator Constructed!")
    
    # Set z to (z + z0)
    def translation(self, z0):
        if type(z0) == type("string"):
            z0 = sympy.sympify(z0, evaluate=False, rational=True, convert_xor=True)
        return PFO(self.dform.subs([(z, z + z0)]), pr=self.pr)
        
    # Set z to (1 / z)
    def translation_inf(self):
        return PFO(to_inv(self.dform, self.deg), pr=self.pr)
    
    # Try to make a list of the special points
    def calclocalexp(self):
        outlst = []
        s = sympy.solve(self.discriminant, z, cubics=True, quartics=True, quintics=True)
        for root in s:
            op = self.translation(root)
            outlst.append((root, op, sympy.roots(op.localind, t)))
        outlst.append((0, self, sympy.roots(self.localind, t)))
        op = self.translation_inf()
        outlst.append((sympy.zoo, op, sympy.roots(op.localind, t)))
        self.localexp = outlst
        return outlst
    
    def __str__(self):
        return "degree = " + str(self.deg) + "\n" + \
               "d-form = " + str(self.dform).replace("**", "^") + "\n" + \
               "t-form = " + str(self.tform).replace("**", "^") + "\n" + \
               "local index = " + str(self.localind).replace("**", "^") + "\n" + \
               "discriminant = " + str(self.discriminant).replace("**", "^")
    
    # Now we are ready to solve the PFeqn. Holomorphic solutions first.
    # Then it's time to find all solutions at the MUM point.
    # Also we want to compute the transition matrices between different points numerically.
    
    # Get the coeff matrix of primtform
    def primtform_list(self):
        if self.primtlist == None:
            outlst = []
            for num in range(self.deg + 1):
                c = self.primtform.coeff(t, num)
                num2 = Z_DEG_BOUND
                while c.coeff(z, num2) == 0:
                    num2 -= 1
                outlst.append([c.coeff(z, num3) for num3 in range(num2 + 1)])
            self.primtlist = outlst
            return outlst
        else:
            return self.primtlist
    
    def hol_sol(self, inlst, termnum):
        plist = self.primtform_list()
        c = [0 for num in range(termnum)]
        outlst = []
        def mono_opr(deg, coeff):
            for num1 in range(self.deg + 1):
                cdn = coeff * (deg ** num1)
                if cdn != 0:
                    for num2 in range(len(plist[num1])):
                        if deg + num2 < termnum:
                            c[deg + num2] += plist[num1][num2] * cdn
                        else:
                            break
        for num in range(termnum):
            den = 0
            for num1 in range(self.deg + 1):
                den += plist[num1][0] * (num ** num1)
            # print(den, c, outlst, den == 0, num)
            if den == 0:
                if num < len(inlst):
                    outlst.append(inlst[num])
                    mono_opr(num, inlst[num])
                else:
                    raise Exception("div 0")
            else:
                if type(den) == type(1):
                    newterm = -c[num] * sympy.Rational(1, den)
                else:
                    newterm = -c[num] / den
                outlst.append(newterm)
                mono_opr(num, newterm)
        return outlst
    
    def log_sol(self, inlst, termnum):
        plist = self.primtform_list()
        c = [0 for num in range(termnum)]
        outlst = []
        def mono_opr(deg, coeff):
            for num1 in range(self.deg + 1):
                cdn = coeff * (deg ** num1)
                if cdn != 0:
                    for num2 in range(len(plist[num1])):
                        if deg + num2 < termnum:
                            c[deg + num2] += plist[num1][num2] * cdn
                        else:
                            break
        def log_mono_opr(deg, logdeg, coeff):
            for num1 in range(self.deg + 1):
                if (deg != 0) or (num1 >= logdeg):
                    cdn = coeff * (deg ** (num1 - logdeg))
                    for num2 in range(num1, num1-logdeg, -1):
                        cdn *= num2
                else:
                    cdn = 0
                if cdn != 0:
                    for num2 in range(len(plist[num1])):
                        if deg + num2 < termnum:
                            c[deg + num2] += cdn * plist[num1][num2]
                        else:
                            break
        for num2 in range(1, len(inlst)):
            for num in range(len(inlst[num2])):
                if num < termnum:
                    log_mono_opr(num, num2, inlst[num2][num])
                else:
                    break
        for num in range(termnum):
            den = 0
            for num1 in range(self.deg + 1):
                den += plist[num1][0] * (num ** num1)
            if den == 0:
                if num < len(inlst[0]):
                    outlst.append(inlst[0][num])
                    mono_opr(num, inlst[0][num])
                else:
                    raise Exception("div 0")
            else:
                if type(den) == type(1):
                    newterm = -c[num] * sympy.Rational(1, den)
                else:
                    newterm = -c[num] / den
                outlst.append(newterm)
                mono_opr(num, newterm)
        return outlst
    
    def all_sol(self, termnum):
        soldict = dict()
        rootsdict = sympy.roots(self.localind, t)
        roots = rootsdict.keys()
        UPPERBOUND = max(roots)
        def listdiv(inlst, num):
            fct = sympy.factorial(num)
            return [term / fct for term in inlst]
        for root in roots:
            # Get the Holomorphic Solution
            initval = [0 for num in range(UPPERBOUND + 1)]
            initval[root] = 1
            sollist = [self.hol_sol(initval, termnum)]
            for logdeg in range(1, rootsdict[root]):
                initval = [[0 for num in range(UPPERBOUND + 1)]]
                for num in range(1, len(sollist) + 1):
                    initval.append(listdiv(sollist[-num], num))
                sollist.append(self.log_sol(initval, termnum))
            soldict[root] = sollist
        return soldict
    
    def isMUM(self):
        return sympy.expand(self.localind - t ** self.deg) == 0
    
    def qcoord(self, termnum):
        def listdiv(inlst, num):
            fct = sympy.Integer(num)
            return [term / fct for term in inlst]
        def mulcoeff(inlst1, inlst2):
            outlst = [0 for num in range(len(inlst1) - 1)]
            for num1 in range(len(inlst1) - 1):
                for num2 in range(len(inlst1) - 1 - num1):
                    outlst[num1 + num2] += inlst1[num1] * inlst2[num2]
            return outlst
        if not self.isMUM():
            raise Exception("Not a local MUM point")
        else:
            holsol = self.hol_sol([1], termnum - 1)
            logsol = self.log_sol([[0], holsol], termnum - 1)[1:]
            mult = logsol
            outstr = [1] + [0 for num in range(termnum - 2)]
            for num in range(1, termnum - 1):
                for num2 in range(num, termnum - 1):
                    outstr[num2] += mult[num2 - num]
                print(mult, outstr)
                mult = listdiv(mult, num + 1)
                mult = mulcoeff(mult, logsol)
            return [0] + outstr
        
if __name__ == "__main__":
    op = PFO(TEST_PFO)
    print(op)
    print(op.all_sol(6))
    
    print()
    
    op3=op.calclocalexp()[0][1]
    print(op3)
    print(op3.all_sol(6))
    
    