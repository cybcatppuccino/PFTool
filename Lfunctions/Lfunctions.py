from LfunctionsL2L4get import L2 as L2
from LfunctionsL2L4get import L4 as L4

THE_BOUND = 10000

def continued_fraction(r, error=1e-7, d=20):
    r1 = r - int(r)
    if d == 0:
        return []
    else:
        if r1 < error:
            return [int(r)]
        else:
            return [int(r)] + continued_fraction(1/r1, error, d-1)

def from_continued_fraction(inlst):
    if len(inlst) == 1:
        return inlst[0], 1
    else:
        rs = from_continued_fraction(inlst[1:])
        return inlst[0] * rs[0] + rs[1], rs[0]

def isrational(r, bound=THE_BOUND, error=1e-7, d=20):
    if r < 0:
        tp = isrational(-r, bound)
        return (-tp[0], tp[1])
    if r == 0:
        return (0, 1)
    else:
        rs = from_continued_fraction(continued_fraction(r, error, d))
        if abs(rs[0]) + abs(rs[1]) > bound:
            return 0, 0
        else:
            return rs
        
def testL21(inval, bound=THE_BOUND, error=1e-7, d=20):
    if inval == 0:
        return []
    else:
        outlst = []
        for Lf in L2:
            if abs(Lf[3]) > 1e-12:
                isr = isrational(inval / Lf[3], bound, error, d)
                if isr[1] != 0:
                    outlst.append((Lf, isr))
        return outlst

def testL41(inval, bound=THE_BOUND, error=1e-7, d=20):
    if inval == 0:
        return []
    else:
        outlst = []
        for Lf in L4:
            if abs(Lf[3]) > 1e-12:
                isr = isrational(inval / Lf[3], bound, error, d)
                if isr[1] != 0:
                    outlst.append((Lf, isr))
        return outlst
    
def testL42(inval, bound=THE_BOUND, error=1e-7, d=20):
    if inval == 0:
        return []
    else:
        outlst = []
        for Lf in L4:
            if abs(Lf[4]) > 1e-12:
                isr = isrational(inval / Lf[4], bound, error, d)
                if isr[1] != 0:
                    outlst.append((Lf, isr))
        return outlst

def alltest(inval, bound=THE_BOUND, error=1e-7, d=20):
    return (testL21(inval, bound, error, d), testL41(inval, bound, error, d), testL42(inval, bound, error, d))
def alltest2(inval, bound=THE_BOUND, error=1e-7, d=20):
    pi = 3.14159265358979323846264338
    zeta3 = 1.20205690315959428539973
    return alltest(inval, bound, error, d)+alltest(inval*pi, bound, error, d)+alltest(inval*(pi**2), bound, error, d)+alltest(inval*(pi**3), bound, error, d)+alltest(inval*zeta3, bound, error, d)

if __name__ == '__main__':
    tl = [1.4325537363690496706,
    0.81107556076091088463,
    -0.82150213771575534351,
    0.19370796305064679945]

    print(alltest2(tl[0]))
    print(alltest2(tl[1]))
    print(alltest2(tl[2]))
    print(alltest2(tl[3]))
    print(alltest2(1.4780203895125552419))
    print(alltest2(1.6180572375999716608))
    print(alltest2(-1.1868969111290801888))
    print(alltest2(0.24591521177202073426))