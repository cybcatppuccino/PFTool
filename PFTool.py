'''
Picard-Fuchs Operator Tool by cybcat

We only deal with PF Operator in polynomial of z, d = d/dz and t = z * d/dz.
'''

import sympy
import fractions

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

DEG_BOUND = 30
Z_DEG_BOUND = 150
# TEST_PFO = -3125*d**4*z**5 + d**4*z**4 - 25000*d**3*z**4 + 6*d**3*z**3 - 45000*d**2*z**3 + 7*d**2*z**2 - 15000*d*z**2 + d*z - 120*z
TEST_PFO = "-3125*d**4*z**5 + d**4*z**4 - 25000*d**3*z**4 + 6*d**3*z**3 - 45000*d**2*z**3 + 7*d**2*z**2 - 15000*d*z**2 + d*z - 120*z"
LEG = "z * (z-1) * d^2 + (2*z - 1) * d + 1/4"


# We introduce
X = sympy.Symbol('t')
# Because on AESZ List they have (e.g. AESZ22): 
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
    while (ineqn.coeff(t, t_deg) == 0) and (t_deg > -1):
        t_deg -= 1
    lst = t_power_to_d_list(t_deg)
    outeqn = 0
    for deg in range(t_deg + 1):
        outeqn += ineqn.coeff(t, deg) * lst[deg]
    return sympy.expand(outeqn)

def to_t_form(ineqn):
    ineqn = sympy.expand(ineqn)
    d_deg = DEG_BOUND
    while (ineqn.coeff(d, d_deg) == 0) and (d_deg > -1):
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
            eqn = sympy.sympify(ineqn.replace("^", "**"), evaluate=False, rational=True, convert_xor=False)
        elif type(ineqn) == type(["List", 114514]):
            eqn = 0
            for num in range(len(ineqn)):
                eqn += ineqn[num] * (z ** num)
        ineqn = sympy.simplify(eqn)
        
        # Standard Forms
        self.dform = to_d_form(ineqn)
        self.tform = to_t_form(ineqn)
        
        # Degree
        self.deg = DEG_BOUND
        while (self.dform.coeff(d, self.deg) == 0) and (self.deg > -1):
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
    
    # [q(x)=e^(log_sol/hol_sol), x(q)]
    def qcoord(self, termnum, pr=False):
        def listdiv(inlst, num):
            fct = sympy.Integer(num)
            return [term / fct for term in inlst]
        def mulcoeff(inlst1, inlst2):
            outlst = [0 for num in range(len(inlst1) - 1)]
            for num1 in range(len(inlst1) - 1):
                for num2 in range(len(inlst1) - 1 - num1):
                    outlst[num1 + num2] += inlst1[num1] * inlst2[num2]
            return outlst
        def mulcoeff2(inlst1, inlst2):
            outlst = [0 for num in range(len(inlst1))]
            for num1 in range(len(inlst1)):
                for num2 in range(len(inlst1) - num1):
                    outlst[num1 + num2] += inlst1[num1] * inlst2[num2]
            return outlst
        if sympy.roots(self.localind, t)[0] < 2:
            raise Exception("t^2 not a factor of local index")
        else:
            if pr:
                print("Start to compute q-coordinates!")
            holsol = self.hol_sol([1], termnum - 1)
            logsol = self.log_sol([[0], holsol], termnum - 1)[1:]
            
            holsol = listdiv(holsol[1:], -1)
            mult = logsol
            for num in range(1, termnum - 1):
                mult = mulcoeff(mult, holsol)
                for num2 in range(num, termnum - 2):
                    logsol[num2] += mult[num2 - num]
            
            mult = logsol
            outlst = [1] + [0 for num in range(termnum - 2)]
            for num in range(1, termnum - 1):
                for num2 in range(num, termnum - 1):
                    outlst[num2] += mult[num2 - num]
                mult = listdiv(mult, num + 1)
                mult = mulcoeff(mult, logsol)
                
            if pr:
                print("The expansion of q(z) in z is done!")
            
            negoutlst = listdiv(outlst[1:], -1)
            mult = negoutlst
            invoutlst = [1] + mult
            for num in range(2, termnum):
                mult = mulcoeff(mult, negoutlst)
                for num2 in range(num, termnum - 1):
                    invoutlst[num2] += mult[num2 - num]
            
            mult = invoutlst
            outlst2 = [1]
            for num in range(2, termnum):
                mult = mulcoeff2(mult, invoutlst)
                outlst2.append(mult[num - 1] / sympy.Integer(num))
            
            if pr:
                print("The expansion of z(q) in q is done!")

            return ([0] + outlst, [0] + outlst2)
    
    # Only CY eqn in the database satisfies some good conditions
    def yukawa(self, termnum, pr=False):
            
        zq = self.qcoord(termnum, pr)[1][1:]
        
        if pr:
            print("Computation of q-coordinates completed!")
        
        yuk = PFO(2 * self.primdform.coeff(d, 4) * d - self.primdform.coeff(d, 3))
        alist = yuk.hol_sol([0, 0, 0, 1], termnum + 3)[3:]
        
        if pr:
            print("ODE a'(z) = g(z)a(z) Solved!")
            print("g(z) = ", str(sympy.simplify(self.primdform.coeff(d,3) / (2 * self.primdform.coeff(d,4)))).replace('**', '^') )
            print("a(z)/z^3 = ", alist)
        
        def listdiv(inlst, num):
            fct = sympy.Integer(num)
            return [term / fct for term in inlst]
        def mulcoeff(inlst1, inlst2):
            outlst = [0 for num in range(len(inlst1) - 1)]
            for num1 in range(len(inlst1) - 1):
                for num2 in range(len(inlst1) - 1 - num1):
                    outlst[num1 + num2] += inlst1[num1] * inlst2[num2]
            return outlst
        def mulcoeff2(inlst1, inlst2):
            outlst = [0 for num in range(len(inlst1))]
            for num1 in range(len(inlst1)):
                for num2 in range(len(inlst2)):
                    if num1 + num2 >= len(inlst1):
                        break
                    else:
                        outlst[num1 + num2] += inlst1[num1] * inlst2[num2]
            return outlst
        
        holsol = self.hol_sol([1], termnum)
        logsol = self.log_sol([[0], holsol], termnum)[1:]
        
        holsol = listdiv(holsol[1:], -1)
        mult = logsol
        for num in range(1, termnum):
            mult = mulcoeff(mult, holsol)
            for num2 in range(num, termnum - 1):
                logsol[num2] += mult[num2 - num]
        
        for num in range(1, termnum):
            logsol[num-1] *= num
        logsol = [1] + logsol
        holsol = [1] + listdiv(holsol, -1)
        
        if pr:
            print("Ready for denominator multiplication!")
        
        mult = mulcoeff2(mulcoeff2(mulcoeff2(mulcoeff2(mulcoeff2(logsol, logsol), logsol), holsol), holsol), alist)
        invlst = [1] + [0 for num in range(termnum - 1)]
        
        zqmult = zq
        for num in range(1, termnum):
            for num2 in range(num, termnum):
                invlst[num2] += mult[num] * zqmult[num2 - num]
            zqmult = mulcoeff(zqmult, zq)
            
        if pr:
            print("Now in q-coordinate!")
            
        outlst = [1] + [0 for num in range(termnum - 1)]
        invlst = listdiv(invlst[1:], -1)
        mult = invlst
        for num in range(1, termnum):
            for num2 in range(num, termnum):
                outlst[num2] += mult[num2 - num]
            mult = mulcoeff(mult, invlst)
        
        return outlst
        
    def instanton(self, termnum, pr=False):
        
        yk = self.yukawa(termnum, pr)
        
        if pr:
            print("Start to compute the instanton!")
        
        outlst = [1]
        for num in range(1, termnum):
            outlst.append(yk[num] / (sympy.Integer(num) ** 3))
            for num2 in range(2, termnum):
                if num2 * num >= termnum:
                    break
                else:
                    yk[num2 * num] -= yk[num]
        
        return outlst
        
if __name__ == "__main__":
    opr = PFO(TEST_PFO)
    print(opr)
    print(opr.all_sol(5))
    print(opr.instanton(5, pr=True))