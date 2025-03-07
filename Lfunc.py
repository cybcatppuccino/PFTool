from LfunctionsL2L4get import L2 as L2
from LfunctionsL2L4get import L4 as L4
from numba import njit
from mpmath import pslq
from math import log

THE_BOUND = 2000
THE_BOUND2 = 100
THE_BOUND3 = 13

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

def isalgebraic(r, bound=THE_BOUND2, error=1e-12):
    theN = round(1/error)
    return pslq([1,r,r*r],tol=error,maxcoeff=bound)

log2 = 0.69314718055994530942
log3 = 1.0986122886681096914
log5 = 1.6094379124341003746
logpi = 1.1447298858494001741

def islogrel(r1, r2, bound=THE_BOUND3, error=1e-12):
    return pslq([log(abs(r1)), log(abs(r2)), log2, log3, log5, logpi], tol=error, maxcoeff=bound)

def testL21L(inval, bound=THE_BOUND3, error=1e-12):
    if abs(inval) < error:
        return []
    else:
        outlst = []
        for Lf in L2[:200]:
            if abs(Lf[3]) > error:
                isr = islogrel(inval, Lf[3], bound, error)
                if isr is not None:
                    if abs(isr[0]) == abs(isr[1]) > 0:
                        outlst.append((Lf, isr))
        return outlst
    
def testL41L(inval, bound=THE_BOUND3, error=1e-12):
    if abs(inval) < error:
        return []
    else:
        outlst = []
        for Lf in L4[:500]:
            if abs(Lf[3]) > error:
                isr = islogrel(inval, Lf[3], bound, error)
                if isr is not None:
                    if abs(isr[0]) == abs(isr[1]) > 0:
                        outlst.append((Lf, isr))
        return outlst
    
def testL42L(inval, bound=THE_BOUND3, error=1e-12):
    if abs(inval) < error:
        return []
    else:
        outlst = []
        for Lf in L4[:500]:
            if abs(Lf[4]) > error:
                isr = islogrel(inval, Lf[4], bound, error)
                if isr is not None:
                    if abs(isr[0]) == abs(isr[1]) > 0:
                        outlst.append((Lf, isr))
        return outlst

def testL21A(inval, bound=THE_BOUND2, error=1e-12):
    if abs(inval) < error:
        return []
    else:
        outlst = []
        for Lf in L2[:200]:
            if abs(Lf[3]) > error:
                isr = isalgebraic(inval / Lf[3], bound, error)
                if isr is not None:
                    outlst.append((Lf, isr))
        return outlst
    
def testL41A(inval, bound=THE_BOUND2, error=1e-12):
    if inval == 0:
        return []
    else:
        outlst = []
        for Lf in L4[:400]:
            if abs(Lf[3]) > 1e-12:
                isr = isalgebraic(inval / Lf[3], bound, error)
                if isr is not None:
                    outlst.append((Lf, isr))
        return outlst
    
def testL42A(inval, bound=THE_BOUND2, error=1e-12):
    if inval == 0:
        return []
    else:
        outlst = []
        for Lf in L4[:400]:
            if abs(Lf[4]) > 1e-12:
                isr = isalgebraic(inval / Lf[4], bound, error)
                if isr is not None:
                    outlst.append((Lf, isr))
        return outlst

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

def alltest2L(inval, bound=THE_BOUND3, error=1e-12):
    return (testL21L(inval, bound, error), testL41L(inval, bound, error), testL42L(inval, bound, error))
def alltestA(inval, bound=THE_BOUND2, error=1e-12):
    return (testL21A(inval, bound, error), testL41A(inval, bound, error), testL42A(inval, bound, error))
def alltest2A(inval, bound=THE_BOUND2, error=1e-12):
    pi = 3.14159265358979323846264338
    return   alltestA(inval, bound, error)\
            +alltestA(inval*pi, bound, error)\
            +alltestA(inval*(pi**2), bound, error)\
            +alltestA(inval/pi, bound, error)\
            +alltestA(inval/(pi**2), bound, error)\

def testlistA(inlst, bound=THE_BOUND2, error=1e-12):
    return [alltest2A(num, bound, error) for num in inlst]

def alltest(inval, bound=THE_BOUND, error=1e-7, d=20):
    return (testL21(inval, bound, error, d), testL41(inval, bound, error, d), testL42(inval, bound, error, d))
def alltest2(inval, bound=THE_BOUND, error=1e-7, d=20):
    pi = 3.14159265358979323846264338
    return   alltest(inval, bound, error, d)\
            +alltest(inval*pi, bound, error, d)\
            +alltest(inval*(pi**2), bound, error, d)\
            +alltest(inval/pi, bound, error, d)\
            +alltest(inval/(pi**2), bound, error, d)\
            +alltest(1/inval, bound, error, d)\
            +alltest(pi/inval, bound, error, d)\
            +alltest((pi**2)/inval, bound, error, d)\
            +alltest(1/(pi*inval), bound, error, d)\
            +alltest(1/((pi**2)*inval), bound, error, d)

def testlist(inlst, bound=THE_BOUND, error=1e-7, d=20):
    return [alltest2(num, bound, error, d) for num in inlst]

def testL41_normal(inval):
    return (testL41(inval), testL41A(inval), testL41L(inval))

def is_alltest2_nontriv(inval, bound=THE_BOUND, error=1e-7, d=20):
    return any(alltest2(inval, bound, error, d))

def getallalgebraicratio():
    pi = 3.14159265358979323846264338
    for Lf in L4:
        l1 = Lf[3]
        l2 = Lf[4]
        if abs(Lf[3]) > 1e-12 and abs(Lf[4]) > 1e-12:
            isr = isalgebraic(l2 / (pi*l1), 20000, 1e-12)
            if isr != None:
                print((Lf, isr))

if __name__ == '__main__':
    #tl = [0.4222831178979863328330828029729991354769372956766793080733248088, - 0.2795357516460809399288839299017105376032984973733782962200]
    #tl = [1.418002734923858,1.99461524031,3.5142998559532]
    #print(testlist(tl, 1000, 1e-8, 15))
    #tl = [-5.9838457209481,-8.761276809822]
    #print(testlist(tl, 1000, 1e-8, 15))
    #tl = [0.6554976628105, 0.664871746772020524,0.3545006837309647, 0.22162391559, 0.41899107750]
    #print(testlist(tl, 1000, 1e-8, 15))
    #getallalgebraicratio()
    for _ in range(5):
        print("-")
    print(testL41_normal(0.78321232657974111623))
    '''
([9, 1, [0, 1, 0, 0, -8, 0, 0, 20, 0, 0, 0, 0, 0], 0.41899107750220405, 0.7599657499933077], [1, 0, -3])
([32, 2, [0, 1, 0, 0, 0, 22, 0, 0, 0, -27, 0, 0, 0], 1.826530444250898, 1.4345536563041807], [1, -8, 16])
([108, 3, [0, 1, 0, 0, 0, 0, 0, 17, 0, 0, 0, 0, 0], 5.985961894730593, 1.8095556464987674], [-1, 0, 108])
([121, 1, [0, 1, 0, 8, -8, 18, 0, 0, 0, 37, 0, 0, -64], 8.416708642567269, 2.65750792952591], [1, 0, -99])
([225, 3, [0, 1, 0, 0, -8, 0, 0, -20, 0, 0, 0, 0, 0], 10.11843933761676, 1.2235212558407853], [1, 0, -675])
([243, 4, [0, 1, 0, 0, -8, 0, 0, -37, 0, 0, 0, 0, 0], 10.458438747785056, 1.0538616417563913], [1, 0, -972])
([256, 1, [0, 1, 0, 10, 0, 0, 0, 0, 0, 73, 0, -18, 0], 19.432481797712505, 3.597343409380413], [1, 0, -288])
([256, 6, [0, 1, 0, 0, 0, -4, 0, 0, 0, -27, 0, 0, 0], 13.032738000433495, 1.7059814149300998], [0, 1, -24])
([288, 1, [0, 1, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0], 15.701877559618492, 2.055370966202662], [0, 1, -24])
([288, 3, [0, 1, 0, 0, 0, -22, 0, 0, 0, 0, 0, 0, 0], 12.654574124055634, 1.1043199195122408], [0, 1, -36])
([432, 13, [0, 1, 0, 0, 0, 0, 0, 37, 0, 0, 0, 0, 0], 24.687237667891367, 2.4876497770400223], [1, 0, -972])
([441, 1, [0, 1, 5, 0, 17, 0, 0, 0, 45, 0, 0, 68, 0], 46.97248961475386, 6.197291470349258], [1, 0, -567])
([576, 1, [0, 1, 0, 0, 0, 22, 0, 0, 0, 0, 0, 0, 0], 35.79254070459021, 3.1234884147459514], [1, -196, 5760])
([576, 14, [0, 1, 0, 0, 0, 0, 0, 20, 0, 0, 0, 0, 0], 31.997314673114413, 2.4181962088040474], [-1, 0, 1728])
([675, 7, [0, 1, 0, 0, -8, 0, 0, 37, 0, 0, 0, 0, 0], 34.604575770463896, 2.092191917735408], [-1, 0, 2700])
([784, 6, [0, 1, 0, 0, 0, 0, 0, 0, 0, -27, 0, 68, 0], 40.67936928537757, 2.0126294429921527], [1, 0, -4032])
([800, 11, [0, 1, 0, 0, 0, 0, 0, 0, 0, -27, 0, 0, 0], 39.08718484278609, 1.5349501593950385], [1, -176, 7680])
([900, 5, [0, 1, 0, 0, 0, 0, 0, 17, 0, 0, 0, 0, 0], 48.74520125796764, 2.357707068031219], [-4, 0, 16875])
([900, 16, [0, 1, 0, 0, 0, 0, 0, -37, 0, 0, 0, 0, 0], 42.00734102510052, 1.269881474074955], [1, 0, -10800])
([972, 3, [0, 1, 0, 0, 0, 0, 0, 20, 0, 0, 0, 0, 0], 54.81772669023711, 2.7618988283199215], [-1, 0, 3888])
    '''