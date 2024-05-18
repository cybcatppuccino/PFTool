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
        if type(ineqn) == type("string"):
            ineqn = sympy.sympify(ineqn, evaluate=False, rational=True, convert_xor=True)
        
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
        self.indpoly = sympy.factor(sympy.cancel(poly / poly.coeff(t, self.deg)))
        
        if pr:
            print("Operator Constructed!")
    
    # Set z to (z + z0)
    def translation(self, z0):
        if type(z0) == type("string"):
            z0 = sympy.sympify(z0, evaluate=False, rational=True, convert_xor=True)
        return PFO(self.dform.subs([(z, z + z0)]), pr=self.pr)
    
    def transfactor(self):
        deltadform = self.dform.subs([(z, z + zdelta)])
        return sympy.factor(to_primitive(to_t_form(deltadform), t, self.deg).subs([(z, 0)]))
        
    # Set z to (1 / z)
    def translation_inf(self):
        pass
    
    def __str__(self):
        return "degree = " + str(self.deg) + "\n" + \
               "dform = " + str(self.dform) + "\n" + \
               "tform = " + str(self.tform) + "\n" + \
               "indpoly = " + str(self.indpoly)
        
if __name__ == "__main__":
    op = PFO(TEST_PFO, pr=True)
    print(op.translation(1))
    print(op.transfactor())
    
    