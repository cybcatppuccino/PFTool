from LfunctionsL2L4get import L2 as L2
from LfunctionsL2L4get import L4 as L4

THE_BOUND = 2000

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
    return   alltest(inval, bound, error, d)\
            +alltest(inval*pi, bound, error, d)\
            +alltest(inval*(pi**2), bound, error, d)\
            +alltest(inval*(pi**3), bound, error, d)\
            +alltest(1/inval, bound, error, d)\
            +alltest(pi/inval, bound, error, d)\
            +alltest((pi**2)/inval, bound, error, d)\
            +alltest((pi**3)/inval, bound, error, d)
def testlist(inlst, bound=THE_BOUND, error=1e-7, d=20):
    return [alltest2(num, bound, error, d) for num in inlst]

def is_alltest2_nontriv(inval, bound=THE_BOUND, error=1e-7, d=20):
    return any(alltest2(inval, bound, error, d))

if __name__ == '__main__':
    tl = [0.894103772609959222407826531437406803285]
    print(testlist(tl))
    print(is_alltest2_nontriv(tl[0]))
