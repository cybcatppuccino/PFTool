import sympy, pickle
from PFTool import *
from time import time

def get_raw_info(num):
    f = open( '.\\info\\' + str(num) + '.txt', 'r', encoding = 'utf-8' )
    outstr = f.read()
    f.close()
    return outstr

TESTSTR = get_raw_info(1)
# Equation List Information: 
def get_info_from_num(num):
    instr = get_raw_info(num)
    outlst = []
    
    num = instr.find('"[') + 1
    num2 = num
    while instr[num2] != ']':
        num2 += 1
    lst = instr[num + 1: num2].replace("X", "t").replace("^", "**").replace("\/", "/").split(',')
    outlst.append([sympy.simplify(term, evaluate=False, rational=True) for term in lst])
    
    num = instr.find('New Number:') + 11
    num2 = num
    while instr[num2] != '&':
        num2 += 1
    outlst.append(instr[num:num2].replace(' ',''))
    
    num = instr.find('AESZ:') + 5
    num2 = num
    while instr[num2] != '&' and instr[num2] != '\\':
        num2 += 1
    outlst.append(instr[num:num2].replace(' ',''))
    
    num = instr.find('Superseeker: <strong>') + 21
    num2 = num
    while instr[num2] != '<':
        num2 += 1
    outlst.append(instr[num:num2].split(' '))
    
    num = instr.find('Note:</h4>') + 10
    num2 = num
    while instr[num2] != '<':
        num2 += 1
    outlst.append(instr[num:num2])
    
    num = instr.find('Yukawa coupling:') + 16
    num2 = num
    while instr[num2] != '<':
        num2 += 1
    outlst.append(instr[num:num2 - 4])
    
    num = instr.find('q-coordinate :') + 14
    num2 = num
    while instr[num2] != '<':
        num2 += 1
    outlst.append(instr[num:num2 - 4])
    
    num = instr.find('q-coordinate :') + 14
    num2 = num
    while instr[num2] != '<':
        num2 += 1
    outlst.append(instr[num:num2 - 4])
    
    return outlst

class AESZ:
    def __init__(self, num):
        info = get_info_from_num(num)
        
        self.num = num
        self.AESZnum = info[2]
        self.newnum = info[1]
        self.ss = (sympy.simplify(info[3][0], evaluate=False, rational=True), sympy.simplify(info[3][1], evaluate=False, rational=True))
        self.note = info[4]
        self.yukawa = [sympy.simplify(term, evaluate=False, rational=True) for term in info[5].split(',')]
        self.qcoord = [sympy.simplify(term, evaluate=False, rational=True) for term in info[6].split(',')]
        
        self.pfo = PFO(info[0])
        # self.pfo.calclocalexp()
    
    def __str__(self):
        return  "number = " + str(self.num) + "\n" + \
                "AESZ number = " + self.AESZnum + "\n" + \
                "new num = " + self.newnum + "\n" + \
                "superseeker = " + str(self.ss) + "\n" + \
                "note = " + self.note + "\n" + \
                "Yukawa coupling = " + str(self.yukawa) + "\n" + \
                "q-coordinate = " + str(self.qcoord) + "\n" + \
                str(self.pfo)

AESZList = []
def get_all_AESZ(a):
    for num in range(1, a+1):
        AESZList.append(AESZ(num))
        print(num, ' finished!')
    with open('AESZpickle.pkl', "wb") as f:  
        pickle.dump(AESZList, f)
        
def load_all_AESZ():
    f = open("AESZpickle.pkl", "rb")
    return pickle.load(f)

def print_all_AESZ():
    outstr = ""
    AESZList = load_all_AESZ()
    for term in AESZList:
        outstr += str(term) + "\n\n"
    f = open("AllAESZinfo.txt", 'w')
    f.write(outstr)
    f.close()
        
if __name__ == '__main__':
    
    # a = AESZ(465)
    a = AESZ(46)
    ap = a.pfo
    print(a.pfo.translation_inf().localind)
    time1 = time()
    print("Start")
    s = ap.all_sol(300)[0]
    time2 = time()
    print(time2 - time1)
    sp = ap.all_sol(150, 23, 15)[0]
    st = ''
    st2 = ''
    for _ in range(150):
        st += str(padic.to_p_adic(s[3][_], 23, 0)[1])
        if type(sp[0][_]) == type(padic.PN(0,0)):
            st2 += str(hex(sp[0][_].acc)[-1])
    print(st)
    print(st2)
    
    '''
    print(ap.all_sol(8))
    print(ap.qcoord(8))
    print(sympy.simplify(ap.primdform.coeff(d,3) / ap.primdform.coeff(d,4)))
    
    Y = sympy.Function('Y')
    C1 = sympy.Symbol('C1')
    ode = 2 * ap.primdform.coeff(d, 4) * Y(z).diff(z) - ap.primdform.coeff(d, 3) * Y(z)
    rst = sympy.cancel(sympy.dsolve(ode).rhs.subs([(C1, 1)]) / (z ** 3))
    ps = sympy.fps(rst).truncate(8).removeO()
    
    print(ps)
    '''

    '''
    b = AESZ(175) # AESZ 34
    c = b.pfo.translation(sympy.Rational(-1, 7))
    print(c, c.all_sol(8))
    '''