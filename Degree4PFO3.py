import AESZ
import PFTool
import sympy
import sys
sys.set_int_max_str_digits(100000)

'''
(2, '', '4.31', 'This is operator "4.31" from ...') (64*z - 1)**2*(108*z - 1)**2
(3, '60', '4.26', 'Sporadic YY-Operator') (16*z - 1)**2*(108*z - 1)**2
(22, '~37', '4.13', 'YY-Operator equivalent to AESZ 37$=C \\ast \\alpha ~tilde c \\ast i$') (256*z - 1)**2*(1024*z - 1)**2
(29, '356\n', '4.32', 'Sporadic YY-Operator.') (108*z - 1)**2*(128*z - 1)**2
(44, '276', '4.55', 'Sporadic Operator.') (110592*z - 1)**2*(331776*z - 1)**2
(74, '362', '4.73', 'Sporadic Operator. There is a second MUM-point corresponding to Operator AESZ 361/4.72') (6912*z - 1)**2*(5308416*z**2 - 3584*z + 1)
(86, '189', '4.27', 'Sporadic YY-Operator') (4*z - 1)**2*(256*z - 1)**2
(90, '~39', '4.1', 'YY-Operator equivalent to AESZ 39=$A \\ast \\alpha$.') (64*z - 1)**2*(256*z - 1)**2
(105, '353', '4.71', 'Sporadic Operator, reducible to 3.33, so not a primary operator.') (16*z - 1)*(64*z - 1)**3
(128, '', '4.7', 'YY-operator equivalent to (:AESZ 50), $\\tilde B \\ast \\alpha$') (108*z - 1)**2*(432*z - 1)**2
(129, '254', '4.49', 'Sporadic Operator.') (4096*z - 1)*(6912*z - 1)*(15360*z - 1)**2
(155, '294', '4.63', 'Sporadic Operator.') (139264*z - 1)**2*(1073741824*z**2 - 22272*z + 1)
(159, '~66', '4.19', 'YY-Operator equivalent to $AESZ 66 =$D \\ast \\alpha \\tilde c \\ast j$') (1728*z - 1)**2*(6912*z - 1)**2
(177, '', '4.77', 'Sporadic Operator.') (16*z - 1)**2*(24*z - 1)*(25*z - 1)
(523, '233', '4.45', 'Sporadic Operator.') (192*z - 1)**2*(432*z - 1)*(512*z - 1)
[2, 3, 22, 29, 44, 74, 86, 90, 105, 128, 129, 155, 159, 177, 523]
'''

TheList = [2, 3, 5, 13, 22, 24, 28, 29, 31, 33, 34, 39, 42, 
           44, 51, 52, 56, 57, 65, 67, 74, 77, 86, 90, 93, 
           94, 96, 99, 105, 116, 117, 118, 121, 128, 129, 130, 
           136, 138, 143, 145, 146, 155, 157, 159, 177, 182, 193, 
           236, 242, 243, 260, 289, 297, 299, 325, 378, 392, 393, 
           398, 399, 420, 421, 432, 439, 441, 445, 470, 478, 488, 
           498, 510, 523, 529, 545, 551]

def get_transed_pfo(innum):
    X = AESZ.AESZ(innum)
    op = X.pfo
    val = min(abs(x) for x in sympy.solve(op.discriminant))
    return op, val

# print(get_transed_pfo(7)[0])

outlist1 = []
outlist2 = []
for n in TheList:
    thepfo, inval = get_transed_pfo(n)
    solinfo = thepfo.all_sol(300)[0]
    outlist1.append(solinfo)
    possiblez = []
    for a in range(-14, 1):
        for b in range(-7, 1):
            for c in range(-2, 1):
                for d in range(-1, 1):
                    for e in range(-1, 1):
                        z0 = (sympy.Integer(2) ** a) * (sympy.Integer(3) ** b) * (sympy.Integer(5) ** c) * (sympy.Integer(7) ** d) * (sympy.Integer(11) ** e)
                        if z0 > 1/10000000 and z0 < inval:
                            possiblez.append(z0)
                            possiblez.append(-z0)
    outlist2.append(possiblez)
    print(str(n) + " done! " + str(len(possiblez)))

st1 = str(outlist1)
st1 = st1.replace('[', '{')
st1 = st1.replace(']', '}')

st2 = str(outlist2)
st2 = st2.replace('[', '{')
st2 = st2.replace(']', '}')

with open('soldegree4pfo3.txt', 'w') as f:
    f.write("sollist = " + st1 + ";\n" + "zlist = " + st2 + ";")

