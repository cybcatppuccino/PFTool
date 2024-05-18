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
    def localexp(self):
        outlst = []
        s = sympy.solve(self.discriminant, z, cubics=True, quartics=True, quintics=True)
        for root in s:
            op = self.translation(root)
            outlst.append((root, op, sympy.roots(op.localind, t)))
        outlst.append((0, self, sympy.roots(self.localind, t)))
        op = self.translation_inf()
        outlst.append((sympy.zoo, op, sympy.roots(op.localind, t)))
        return outlst
    
    def __str__(self):
        return "deg = " + str(self.deg) + "\n" + \
               "dform = " + str(self.dform).replace("**", "^") + "\n" + \
               "tform = " + str(self.tform).replace("**", "^") + "\n" + \
               "localind = " + str(self.localind).replace("**", "^") + "\n" + \
               "discriminant = " + str(self.discriminant).replace("**", "^")
        
if __name__ == "__main__":
    # op = PFO(TEST_PFO)
    # print(op.localexp())
    op2 = PFO(TEST_LIST22)
    print(op2)
    print(op2.localexp())
    '''
    op2 = PFO(LEG)
    print(op2)
    print(op2.discriminant())
    '''
    
    