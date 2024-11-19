import AESZ
import sympy

import sys
sys.set_int_max_str_digits(100000)

'''
Attractor Point
level N = 27
infty index a1, a2, a3, a4 = 1/3, 1/3, 2/3, 2/3
conifold pt = 1/3^6 = 1/729
mirror M = X_{3,3}(1^6)

kappa = 9
c2.D = 54
Euler chi(M) = -144

AESZ: (207, '4', '1.4', 'A-incarnation: X(3,3) in P^5.')
'''

X = AESZ.AESZ(207)
op = X.pfo
print(op)
# print(op.calclocalexp())
# (1/729, {2: 1, 1: 2, 0: 1}), (0, {0: 4}), (zoo, {2/3: 2, 1/3: 2})

'''
solinfo = op.all_sol(2000)

with open('solin.txt', 'w') as f:
    f.write(str(solinfo))
'''
    
'''
(23, '393', '3.23', 'This is operator "3.23" from ...')
(46, '', '3.24', 'This is operator $\\tilde{C_9}$')
(47, '', '3.25', 'This is operator $\\tilde{C_17}$')
(89, '', '3.15', 'Operator equivalent to AESZ 328')
(135, '408', '3.27', 'This is operator "3.27" from ...')
(142, '228', '3.3', 'This is operator "3.3" from ...')
(162, '', '3.4', 'Operator equivalent to AESZ 165= $f \\ast f$.')
(165, '227', '3.2', '>console.log("JSON Response:")')
(175, '34', '3.1', '>console.log("JSON Response:")')
(185, '386', '3.16', '>console.log("JSON Response:")')
(188, '~103', '3.10', 'Operator equivalent to AESZ 103 =$c \\ast c$.')
(219, '422', '3.30', 'This is operator "3.30" from ...')
(224, '', '3.32', 'Operator equivalent to AESZ 220')
(227, '~101', '3.9', 'Operator equivalent to $AESZ 101=$b \\ast b$.')
(247, '', '3.34', 'Operator equivalent to AESZ 107 $=d \\ast d$')
(254, '392', '3.22', 'This is operator "3.22" from ...')
(272, '390', '3.20', 'This is operator "3.20" from ...')
(278, '~33', '3.6', 'Operator AESZ 33 is replaced by this equivalent operator.')
(284, '407', '3.26', 'This is operator "3.26" from ...')
(334, '', '3.12', 'Operator equivalent to AESZ 154')
(338, '', '3.5', 'Operator equivalent to AESZ 214')
(342, '388', '3.18', 'This is operator "3.18" from ...')
(359, '411', '3.29', 'This is operator "3.29" from ...')
(384, '', '3.13', 'Operator equivalent to AESZ 229')
(395, '', '3.11', 'Operator equivalent to AESZ 144=c \\ast c$')
(408, '', '3.33', 'Operator equivalent to AESZ 353')
(410, '', '3.31', 'This is operator "3.31" from ...')
(425, '', '3.14', 'This is operator Pi = 3.14 (approx.), equivalent to AESZ 238.')
(447, '~100', '3.8', 'Operator equivalent to AESZ 100= $ a \\ast a$')
(489, '389', '3.19', 'This is operator "3.19" from ...')
(493, '387', '3.17', 'This is operator "3.17" from ...')
(500, '410', '3.28', 'This is operator "3.28" from ...')
(539, '~73', '3.7', 'Operator equivalent to AESZ 73')
(550, '391', '3.21', 'This is operator "3.21" from ...')
'''

# TryList = [23, 46, 47, 89, 135, 142, 162, 165, 175, 185, 188, 219, 224, 227, 247, 254, 272, 278, 284, 334, 338, 342, 359, 384, 395, 408, 425, 447, 489, 493, 500, 539, 550]

TryList = [2, 3, 5, 13, 22, 24, 28, 29, 31, 33, 34, 39, 42, 
           44, 51, 52, 56, 57, 65, 67, 74, 77, 86, 90, 93, 
           94, 96, 99, 105, 116, 117, 118, 121, 128, 129, 130, 
           136, 138, 143, 145, 146, 155, 157, 159, 177, 182, 193, 
           236, 242, 243, 260, 289, 297, 299, 325, 378, 392, 393, 
           398, 399, 420, 421, 432, 439, 441, 445, 470, 478, 488, 
           498, 510, 523, 529, 545, 551]

lenTryList = len(TryList)

outlist1 = []
outlist2 = []
for n in range(lenTryList):
    AESZnum = TryList[n]
    thepfo = AESZ.AESZ(AESZnum).pfo
    thebound = min(abs(_) for _ in sympy.solve(thepfo.discriminant))
    solinfo = thepfo.all_sol(50)[0]
    possiblez = []
    for a in range(-20, 9):
        for b in range(-24, 9):
            for c in range(-6, 4):
                z0 = (sympy.Integer(2) ** a) * (sympy.Integer(3) ** b) * (sympy.Integer(5) ** c)
                if thebound / 100 <= z0 <= thebound / 2:
                    possiblez.append(z0)
    outlist1.append(solinfo)
    outlist2.append(possiblez)
    print(str(n) + " done! " + str(len(possiblez)))

st1 = str(outlist1)
st1 = st1.replace('[', '{')
st1 = st1.replace(']', '}')

st2 = str(outlist2)
st2 = st2.replace('[', '{')
st2 = st2.replace(']', '}')

with open('solin2.txt', 'w') as f:
    f.write("sollist = " + st1 + ";\n" + "zlist = " + st2 + ";")
