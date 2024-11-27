'''
Picard-Fuchs Operator Tool by cybcat

We only deal with PF Operator in polynomial of z, d = d/dz and t = z * d/dz.
'''
import padic
import sympy
from LLL import reduction
import mpmath
mpmath.mp.dps = 64
import numpy as np
import sys
sys.set_int_max_str_digits(100000)

# The Symbols
# z: the variable of the ODE
# d: the differential operator d/dz
# t: the differential operator z * d/dz

z = sympy.Symbol('z')
d = sympy.Symbol('d')
t = sympy.Symbol('t')

DEG_BOUND = 30
Z_DEG_BOUND = 150
TEST_PFO = "-3125*d^4*z^5 + d^4*z^4 - 25000*d^3*z^4 + 6*d^3*z^3 - 45000*d^2*z^3 + 7*d^2*z^2 - 15000*d*z^2 + d*z - 120*z"
TEST_PFO2 = "65536*t**4*z**2 - 512*t**4*z + t**4 + 262144*t**3*z**2 - 1024*t**3*z + 360448*t**2*z**2 - 832*t**2*z + 196608*t*z**2 - 320*t*z + 36864*z**2 - 48*z"
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
    while (not ineqn.coeff(t, t_deg)) and (t_deg > -1):
        t_deg -= 1
    lst = t_power_to_d_list(t_deg)
    outeqn = 0
    for deg in range(t_deg + 1):
        outeqn += ineqn.coeff(t, deg) * lst[deg]
    return sympy.expand(outeqn)

def to_t_form(ineqn):
    ineqn = sympy.expand(ineqn)
    d_deg = DEG_BOUND
    while (not ineqn.coeff(d, d_deg)) and (d_deg > -1):
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

def toc(innum):
    if str(type(innum))[8:14] == "mpmath":
        return mpmath.mpc(innum)
    else:
        repart = sympy.N(sympy.re(innum), mpmath.mp.dps+10)
        impart = sympy.N(sympy.im(innum), mpmath.mp.dps+10)
        return mpmath.mpf(str(repart)) + mpmath.mpf(str(impart)) * mpmath.j

def dlist(inlst):
    return [inlst[num] * num for num in range(1, len(inlst))]

def polyval(inlst, pt):
    return mpmath.polyval(list(reversed(inlst)), toc(pt))

def wrval(inlst, diff, pt):
    outlst = []
    inlstcopy = inlst.copy()
    pt0 = toc(pt)
    for diffnum in range(diff + 1):
        outlst.append(mpmath.polyval(list(reversed(inlstcopy)), pt0))
        if diffnum < diff:
            inlstcopy = dlist(inlstcopy)
    return outlst

def mptomma(inm):
    outstr = '{'
    for num in range(inm.rows):
        outstr += '{'
        for num2 in range(inm.cols):
            outstr += str(inm[num, num2])
            if num2 < inm.cols - 1:
                outstr += ','
        if num < inm.rows - 1:
            outstr += '},'
        else:
            outstr += '}'
    outstr += '}'
    return outstr.replace('j', 'I').replace('.0 ', ' ').replace('.0I', 'I').replace('e', '*10^')

STANDWR = mpmath.matrix([[toc(int(num1 == num2) * sympy.factorial(num1))
                     for num1 in range(4)] for num2 in range(4)])

STANDCONTOUR1 = [0.25000, 0.31250, 0.39063, 0.48120, 0.58496, 0.66797, 0.73438, 0.78750, 0.83000]
STANDCONTOUR2 = [0.83000, 0.83371 - 0.03534 *sympy.I, 0.84470 - 0.06915 *sympy.I, 0.86247 - \
 0.09992 *sympy.I, 0.88625 - 0.12633 *sympy.I, 0.91500 - 0.14722 *sympy.I, 0.94747 - \
 0.16168 *sympy.I, 0.98223 - 0.16907 *sympy.I, 1.01777 - 0.16907 *sympy.I, 1.05253 - \
 0.16168 *sympy.I, 1.08500 - 0.14722 *sympy.I, 1.11375 - 0.12633 *sympy.I, 1.13753 - \
 0.09992 *sympy.I, 1.15530 - 0.06915 *sympy.I, 1.16629 - \
 0.03534 *sympy.I, 1.1700, 1.16629 + 0.03534 *sympy.I, 1.15530 + \
 0.06915 *sympy.I, 1.13753 + 0.09992 *sympy.I, 1.11375 + 0.12633 *sympy.I, 1.08500 + \
 0.14722 *sympy.I, 1.05253 + 0.16168 *sympy.I, 1.01777 + 0.16907 *sympy.I, 0.98223 + \
 0.16907 *sympy.I, 0.94747 + 0.16168 *sympy.I, 0.91500 + 0.14722 *sympy.I, 0.88625 + \
 0.12633 *sympy.I, 0.86247 + 0.09992 *sympy.I, 0.84470 + 0.06915 *sympy.I, 0.83371 + \
 0.03534 *sympy.I, 0.83000]

def approx(inval, extreprec=20):
    r = abs(inval)
    if r == 0:
        return 0
    else:
        d = 2 ** (sympy.floor(sympy.log(r, 2)) - extreprec)
        #print(inval, d, sympy.floor(inval / d) * d)
        return sympy.simplify(sympy.floor(inval / d) * d)

def hololist(inpfo, zlist, termnum, pr=False):
    if pr:
        print(inpfo)
        print(zlist)
    outmatrix = mpmath.eye(4)
    for num in range(len(zlist) - 1):
        newpfo = inpfo.translation(zlist[num])
        outmatrix = newpfo.CY_trans_holbasis(termnum, toc(zlist[num + 1]-zlist[num])) * outmatrix
        if pr:
            print(num+1, "/", len(zlist) - 1)
    return outmatrix

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
        while (not self.dform.coeff(d, self.deg)) and (self.deg > -1):
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
        self.mpprimtlist = None
        self.primtlistpadic = None
        self.mpallsol = None
        
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
    def primtform_list(self, inp=None, inacc=4):
        if inp == None:
            if self.primtlist == None:
                outlst = []
                for num in range(self.deg + 1):
                    c = self.primtform.coeff(t, num)
                    num2 = Z_DEG_BOUND
                    while (not c.coeff(z, num2)):
                        num2 -= 1
                    outlst.append([c.coeff(z, num3) for num3 in range(num2 + 1)])
                self.primtlist = (None, None, outlst)
                return (None, None, outlst)
            else:
                return self.primtlist
        else:
            if self.primtlistpadic == None or self.primtlistpadic[0] != inp or self.primtlistpadic[1] != inacc:
                outlst = []
                for num in range(self.deg + 1):
                    c = self.primtform.coeff(t, num)
                    num2 = Z_DEG_BOUND
                    while not c.coeff(z, num2):
                        num2 -= 1
                    outlst.append([padic.PN(inp, inacc).setval(c.coeff(z, num3)) for num3 in range(num2 + 1)])
                self.primtlistpadic = (inp, inacc, outlst)
                return (inp, inacc, outlst)
            else:
                return self.primtlistpadic
    
    def mpprimtform_list(self):
        if self.mpprimtlist == None:
            rs = self.primtform_list()
            outlst = []
            for num in range(len(rs[2])):
                outlst.append([toc(_) for _ in rs[2][num]])
            self.mpprimtlist = outlst
            return outlst
        else:
            return self.mpprimtlist
    
    def hol_sol(self, inlst, termnum, inp=None, inacc=4, expnd=True):
        plist = self.primtform_list(inp, inacc)[2]
        c = [0 for num in range(termnum)]
        outlst = []
        def mono_opr(deg, coeff):
            for num1 in range(self.deg + 1):
                cdn = coeff * (deg ** num1)
                if expnd:
                    cdn = sympy.expand(cdn)
                if cdn:
                    for num2 in range(len(plist[num1])):
                        if deg + num2 < termnum:
                            c[deg + num2] += plist[num1][num2] * cdn
                            if expnd:
                                c[deg + num2] = sympy.expand(c[deg + num2])
                        else:
                            break
        for num in range(termnum):
            den = 0
            for num1 in range(self.deg + 1):
                den += plist[num1][0] * (num ** num1)
                if expnd:
                    den = sympy.expand(den)
            # print(den, c, outlst, not den, num)
            if not den:
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
                if expnd:
                    newterm = sympy.expand(newterm)
                outlst.append(newterm)
                mono_opr(num, newterm)
        return outlst
    
    def mp_hol_sol(self, inlst, termnum):
        plist = self.mpprimtform_list()
        c = [toc(0) for _ in range(termnum)]
        outlst = [toc(0) for _ in range(termnum)]
        inlst = [toc(_) for _ in inlst]
        def mono_opr(deg, coeff):
            for num1 in range(self.deg + 1):
                cdn = coeff * (deg ** num1)
                if cdn:
                    for num2 in range(len(plist[num1])):
                        if deg + num2 < termnum:
                            c[deg + num2] += plist[num1][num2] * cdn
                        else:
                            break
        for num in range(termnum):
            den = mpmath.mpc(0)
            for num1 in range(len(plist)):
                den += plist[num1][0] * (num ** num1)
            if num < len(inlst):
                outlst[num] = inlst[num]
                mono_opr(num, inlst[num])
            else:
                newterm = -c[num] / den
                outlst[num] = newterm
                mono_opr(num, newterm)
        return outlst
    
    def log_sol(self, inlst, termnum, inp=None, inacc=4):
        plist = self.primtform_list(inp, inacc)[2]
        c = [0 for num in range(termnum)]
        outlst = []
        def mono_opr(deg, coeff):
            for num1 in range(self.deg + 1):
                cdn = coeff * (deg ** num1)
                if cdn :
                    for num2 in range(len(plist[num1])):
                        if deg + num2 < termnum:
                            c[deg + num2] += plist[num1][num2] * cdn
                        else:
                            break
        def log_mono_opr(deg, logdeg, coeff):
            for num1 in range(self.deg + 1):
                if (deg) or (num1 >= logdeg):
                    cdn = coeff * (sympy.Integer(deg) ** (num1 - logdeg))
                    for num2 in range(num1, num1-logdeg, -1):
                        cdn *= num2
                else:
                    cdn = 0
                if cdn :
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
            if not den:
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
    
    def mp_log_sol(self, inlst, termnum):
        plist = self.mpprimtform_list()
        c = [toc(0) for _ in range(termnum)]
        outlst = [toc(0) for _ in range(termnum)]
        inlst = [[toc(_) for _ in __] for __ in inlst]
        def mono_opr(deg, coeff):
            for num1 in range(self.deg + 1):
                cdn = coeff * (deg ** num1)
                if cdn :
                    for num2 in range(len(plist[num1])):
                        if deg + num2 < termnum:
                            c[deg + num2] += plist[num1][num2] * cdn
                        else:
                            break
        def log_mono_opr(deg, logdeg, coeff):
            for num1 in range(self.deg + 1):
                if (deg) or (num1 >= logdeg):
                    cdn = coeff * (deg ** (num1 - logdeg))
                    for num2 in range(num1, num1-logdeg, -1):
                        cdn *= num2
                else:
                    cdn = 0
                if cdn :
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
            if num < len(inlst[0]):
                    outlst[num] = inlst[0][num]
                    mono_opr(num, inlst[0][num])
            else:
                newterm = -c[num] / den
                outlst[num] = newterm
                mono_opr(num, newterm)
        return outlst
    
    def all_sol(self, termnum, inp=None, inacc=4):
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
            sollist = [self.hol_sol(initval, termnum, inp, inacc)]
            for logdeg in range(1, rootsdict[root]):
                initval = [[0 for num in range(UPPERBOUND + 1)]]
                for num in range(1, len(sollist) + 1):
                    initval.append(listdiv(sollist[-num], num))
                sollist.append(self.log_sol(initval, termnum, inp, inacc))
            soldict[root] = sollist
        return soldict
    
    def mp_all_sol(self, termnum):
        if self.mpallsol == None or self.mpallsol[0] < termnum:
            soldict = dict()
            rootsdict = sympy.roots(self.localind, t)
            roots = rootsdict.keys()
            UPPERBOUND = max(roots)
            def listdiv(inlst, num):
                fct = sympy.factorial(num)
                return [term / toc(fct) for term in inlst]
            for root in roots:
                # Get the Holomorphic Solution
                initval = [0 for num in range(UPPERBOUND + 1)]
                initval[root] = 1
                sollist = [self.mp_hol_sol(initval, termnum)]
                for logdeg in range(1, rootsdict[root]):
                    initval = [[0 for num in range(UPPERBOUND + 1)]]
                    for num in range(1, len(sollist) + 1):
                        initval.append(listdiv(sollist[-num], num))
                    sollist.append(self.mp_log_sol(initval, termnum))
                soldict[root] = sollist
            self.mpallsol = (termnum, soldict)
            return soldict
        else:
            return self.mpallsol[1]
    
    def isMUM(self):
        return not sympy.expand(self.localind - t ** self.deg)

    def eval_MUM(self, termnum, pt):
        soldict = self.mp_all_sol(termnum)
        logpt = mpmath.log(toc(pt))
        f0 = polyval(soldict[0][0], pt)
        f1 = polyval(soldict[0][1], pt)
        f2 = polyval(soldict[0][2], pt)
        f3 = polyval(soldict[0][3], pt)
        return [f0, f1+logpt*(f0), f2+logpt*(f1+logpt*(f0/2)), f3+logpt*(f2+logpt*(f1/2+logpt*(f0/6)))]

    def attrval_MUM(self, termnum, pt):
        soldict = self.mp_all_sol(termnum)
        logpt = mpmath.log(toc(pt))
        f0 = polyval(soldict[0][0], pt)
        f1 = polyval(soldict[0][1], pt)
        f2 = polyval(soldict[0][2], pt)
        return mpmath.re(logpt ** 2 + (2*f1*logpt + 2*f2) / f0) / (mpmath.pi ** 2)

    def Wronskian_MUM(self, termnum, pt):
        pt = toc(pt)
        soldict = self.mp_all_sol(termnum)
        logpt = mpmath.log(pt)
        f0 = wrval(soldict[0][0], 3, pt)
        f1 = wrval(soldict[0][1], 3, pt)
        f2 = wrval(soldict[0][2], 3, pt)
        f3 = wrval(soldict[0][3], 3, pt)
        r00 = f0[0]
        r01 = f0[1]
        r02 = f0[2]
        r03 = f0[3]
        r10 = logpt*f0[0] + f1[0]
        r11 = f0[0]/pt + logpt*f0[1] + f1[1]
        r12 = -((f0[0] - 2*pt*f0[1])/pt**2) + logpt*f0[2] + f1[2]
        r13 = (2*f0[0] - 3*pt*f0[1] + 3*pt**2*f0[2])/pt**3 + logpt*f0[3] + f1[3]
        r20 = logpt**2*f0[0]/2 + logpt*f1[0] + f2[0]
        r21 = (logpt*f0[0])/pt + logpt**2*f0[1]/2 + f1[0]/pt + logpt*f1[1] + f2[1]
        r22 = (1/(2*pt**2))*(-2*(-1 + logpt)*f0[0] - 2*f1[0] + pt*(4*f1[1] + logpt*(4*f0[1] + logpt*pt*f0[2] + 2*pt*f1[2]) + 2*pt*f2[2]))
        r23 = (1/(2*pt**3))*((-6 + 4*logpt)*f0[0] + 4*f1[0] + pt*(-6*(-1 + logpt)*f0[1] - 6*f1[1] + pt*(6*f1[2] + logpt*(6*f0[2] + logpt*pt*f0[3] + 2*pt*f1[3]) + 2*pt*f2[3])))
        r30 = logpt**2*(logpt*f0[0] + 3*f1[0])/6 + logpt*f2[0] + f3[0]
        r31 = (1/(6*pt))*(6*f2[0] + logpt*(6*f1[0] + logpt*(3*f0[0] + logpt*pt*f0[1] + 3*pt*f1[1]) + 6*pt*f2[1])) + f3[1]
        r32 = (1/(6*pt**2))*(-3*(-2 + logpt)*logpt*f0[0] - 6*(-1 + logpt)*f1[0] - 6*f2[0] + pt*(12*f2[1] + logpt*(12*f1[1] + logpt*(6*f0[1] + logpt*pt*f0[2] + 3*pt*f1[2]) + 6*pt*f2[2]) + 6*pt*f3[2]))
        r33 = (1/(6*pt**3))*(6*(1 + (-3 + logpt)*logpt)*f0[0] + 6*(-3 + 2*logpt)*f1[0] + 12*f2[0] + pt*(-9*(-2 + logpt)*logpt*f0[1] - 18*(-1 + logpt)*f1[1] - 18*f2[1] + pt*(18*f2[2] + logpt*(18*f1[2] + logpt*(9*f0[2] + logpt*pt*f0[3] + 3*pt*f1[3]) + 6*pt*f2[3]) + 6*pt*f3[3])))
        return mpmath.matrix([[r00,r01,r02,r03],[r10,r11,r12,r13],[r20,r21,r22,r23],[r30,r31,r32,r33]])
    
    def Wronskian_General(self, termnum, pt):
        pt = toc(pt)
        soldict = self.mp_all_sol(termnum)
        f0 = wrval(soldict[0][0], 3, pt)
        f1 = wrval(soldict[1][0], 3, pt)
        f2 = wrval(soldict[2][0], 3, pt)
        f3 = wrval(soldict[3][0], 3, pt)
        return mpmath.matrix([f0, f1, f2, f3])
    
    def Wronskian0(self, termnum, pt):
        if self.isMUM():
            return self.Wronskian_MUM(termnum, pt)
        else:
            return self.Wronskian_General(termnum, pt)

    def Wronskian0_MMA(self, termnum, pt):
        wrsk = self.Wronskian0(termnum, pt)
        outstr = '{'
        for num in range(5):
            if num == 0:
                mult = mpmath.zeta(3)
            else:
                mult = (2 * mpmath.pi * mpmath.j) ** (4-num)
            outstr += '{'
            for num2 in range(4):
                outstr += str(wrsk[num-(num>0), num2]*mult)
                if num2 < 3:
                    outstr += ','
            if num < 4:
                outstr += '},'
            else:
                outstr += '}'
        outstr += '}'
        return outstr.replace('j', 'I').replace('.0 ', ' ').replace('.0I', 'I').replace('e', '*10^')
    
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
    
    def CY_trans_holbasis(self, termnum, znew):
        # newbasis [p0', p1', p2', p3'] = A . oldbasis [p0, p1, p2, p3]
        # using Wronskian computations
        
        # znew need to be small to maintain convergence
        return STANDWR * (self.Wronskian0(termnum, znew) ** -1)

    def monodromy_only_pt(self, pt, termnum = round(mpmath.mp.dps * 2), pr=False):
        # Compute the monodromy to the pt, where pt should be the only special point other than the origin and infinity
        # Return the monodromy matrix under basis at the MUM basis
        contour1 = [approx(_ * pt) for _ in STANDCONTOUR1]
        contour2 = [approx(_ * pt) for _ in STANDCONTOUR2]
        if pr:
            print(termnum)
        mat1 = hololist(self, contour1, termnum, pr)
        mat2 = hololist(self, contour2, termnum, pr)
        return (mat1 ** -1) * (mat2 ** -1) * mat1

    def monodromy_origin_test(self, r, termnum = round(mpmath.mp.dps * 2), pr=False):
        # Compute the monodromy at origin
        # Return the monodromy matrix under basis at the MUM basis
        contour1 = [approx(_ * r) for _ in [0, 1]]
        contour2 = [approx(sympy.exp(2*sympy.pi*sympy.I*_/25) * r) for _ in range(26)]
        if pr:
            print(termnum)
        mat1 = hololist(self, contour1, termnum, pr)
        mat2 = hololist(self, contour2, termnum, pr)
        return (mat1 ** -1) * (mat2 ** -1) * mat1

    def CY_trans_holbasis_real_avoidreal(self, pts, znew, termnum = round(mpmath.mp.dps * 2), pr=False):
        contour1 = [0]
        contour2 = [znew]
        pts0 = pts.copy()
        def set_dist(inpts, z0):
            return min(abs(z0-_) for _ in inpts)
        while sympy.im(contour1[-1]) < abs(znew) * 1:
            contour1.append(approx(contour1[-1] + set_dist(pts0, contour1[-1]) * sympy.I / 5))
            if 0 not in pts0:
                pts0.append(0)
        while sympy.im(contour2[-1]) < sympy.im(contour1[-1]):
            contour2.append(approx(contour2[-1] + set_dist(pts0, contour2[-1]) * sympy.I / 5))
        
        end1 = contour1[-1]
        end2 = contour2[-1]
        contour3 = [approx((end1 * (6-_) + end2 * _)/6) for _ in range(1,6)]
        contour = contour1 + contour3 + list(reversed(contour2))
        
        #contour = contour1 + list(reversed(contour2))
        return hololist(self, contour, termnum, pr)

    def eval_MUM0(self, znew, pr=False):
        pts = list(sympy.solve(self.discriminant))
        mat = self.CY_trans_holbasis_real_avoidreal(pts, znew=znew, pr=pr)
        return (mat ** -1) * mpmath.matrix([[1], [0], [0], [0]])
    
    def eval_Wronskian0(self, znew, pr=False):
        pts = list(sympy.solve(self.discriminant))
        return self.CY_trans_holbasis_real_avoidreal(pts, znew=znew, pr=pr)
    
    def eval_attr_LLL(self, znew, pr=False):
        pts = list(sympy.solve(self.discriminant))
        mat = self.CY_trans_holbasis_real_avoidreal(pts, znew=znew, pr=pr) ** -1
        sol1 = list(mat * mpmath.matrix([[1], [0], [0], [0]]))
        sol2 = list(mat * mpmath.matrix([[0], [1], [0], [0]]))

        tpj = 2 * mpmath.pi * mpmath.j
        vals = [mpmath.zeta(3), tpj ** 3, tpj ** 2, tpj ** 1]
        tol = 10**(12 - mpmath.mp.dps)
        ls2 = [mpmath.re(sol1[0]*vals[0]), mpmath.re(sol1[0]*vals[1]), mpmath.re(sol1[1]*vals[2]), mpmath.re(sol1[2]*vals[3]), mpmath.re(sol1[3]), \
               mpmath.im(sol1[0]*vals[0]), mpmath.im(sol1[0]*vals[1]), mpmath.im(sol1[1]*vals[2]), mpmath.im(sol1[2]*vals[3]), mpmath.im(sol1[3])]
        ls3 = []
        for _ in ls2:
            if abs(_) > tol:
                ls3.append(_)
        l = len(ls3)
        llllst = [[int(a==b) + int(ls3[a] * 10**(mpmath.mp.dps - 12)) * int(b==l) for b in range(l+1)] for a in range(l)]
        l = 5
        llllst2 = [[int(a==b) + int(ls2[a] * 10**(mpmath.mp.dps - 12)) * int(b==l) + int(ls2[a + l] * 10**(mpmath.mp.dps - 12)) * int(b==l+1) for b in range(l+2)] for a in range(l)]
        term1 = sum(mpmath.conj(sol1[_]) * sol2[3-_] for _ in range(4))
        term2 = sum(mpmath.conj(sol1[_]) * sol1[3-_] for _ in range(4))
        term3 = mpmath.zeta(3) * mpmath.conj(sol1[0]) * sol2[0]
        term4 = mpmath.zeta(3) * mpmath.conj(sol1[0]) * sol1[0]
        ls2 = [mpmath.re(term1), mpmath.re(term2), mpmath.re(term3), mpmath.re(term4), mpmath.im(term1), mpmath.im(term2), mpmath.im(term3), mpmath.im(term4)]
        l = 4
        llllst3 = [[int(a==b) + int(ls2[a] * 10**(mpmath.mp.dps - 12)) * int(b==l) + int(ls2[a + l] * 10**(mpmath.mp.dps - 12)) * int(b==l+1) for b in range(l+2)] for a in range(l)]
        if pr:
            print("LLL")
        return mat, ls2, reduction(llllst), reduction(llllst2), reduction(llllst3)
        
if __name__ == "__main__":
    '''
    opr = PFO(TEST_PFO2)
    z0 = sympy.Integer(-1)/768
    print("Start!")
    opr.mp_all_sol(1000)
    print("Finish!")
    print(STANDWR)
    print(opr.attrval_MUM(1000, z0))
    print(mptomma(opr.CY_trans_holbasis(1000, z0)))
    '''

    b = PFO(TEST_PFO)

    #print(mptomma(b.monodromy_only_pt(sympy.Integer(1)/3125, termnum=1000, pr=True)))
    #z0 = sympy.Integer(20)/3125
    #print(b.eval_MUM0(z0, pr=True))

    