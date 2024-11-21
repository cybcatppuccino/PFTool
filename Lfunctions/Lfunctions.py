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
    zeta3 = 1.20205690315959428539973
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

if __name__ == '__main__':
    tl = [-11.054835954699885689, 14.992300192559700866, 
        -6.0093536649054504750, 0.91142256061883843065, 
        56.519640174430078546, -18.878941326489863311, 2.8633184207561408108,1.4325537363690496706, 0.81107556076091088463, 
        -0.82150213771575534351, 0.19370796305064679945, 
        4.5501232557012492971, -2.5808250807561275957, 0.60855151366175509889,
        6.6289922379644914060, -1.0725245331139746826, \
        -1.6874079929248544282, 0.65201673108332164257, \
        3.3694351940256857923, -5.3011485541814204776, 2.0483709723889950625
    ]

    print(testlist(tl))
