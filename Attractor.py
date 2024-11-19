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
(7, '143', '2.30', '$D\\ast h^{\\tilde{\\;}}B\\ast\\kappa$')
(8, '~82', '2.45', 'Operator equivalent to $\\hat{8}$')
(17, '406', '2.68', 'This is operator "2.68" from ...')
(27, '~94', '2.48', 'Operator equivalent to $\\hat{11}$')
(38, '', '2.33', 'This is operator "2.33" from ...')
(45, '', '2.66', 'This is operator "2.66" from ...')
(54, '', '2.50', 'Operator equivalent to $\\hat{13}$ of AESZ')
(66, '~80,~81', '2.39', 'Operator equivalent to $\\hat{1}$ of AESZ')
(68, '', '2.37', 'Hadamard product $B\\ast c$.')
(73, '', '2.49', 'Operator equivalent to $\\hat{12}$')
(75, '', '2.46', 'Operator equivalent to $\\hat{9}$')
(81, '~98', '2.42', 'Operator equivalent to $\\hat{5}$')
(82, '112', '2.19', '>console.log("JSON Response:")')
(98, '', '2.47', 'Operator equivalent to $\\hat{10}$')
(100, '', '2.36', 'This is operator "2.36" from ...')
(102, '47', '2.59', 'Hadamard product $I \\ast \\kappa$')
(104, '', '2.34', 'This is operator "2.34" from ...')
(106, '84', '2.63', 'This is operator "2.63" from ...')
(107, '183', '2.65', 'This is operator "2.65" from ...')
(113, '61', '2.38', 'This is operator "2.38" from ...')
(126, '', '2.41', 'Operator equivalent to $\\hat{4}$')
(131, '~88,~89', '2.51', 'Operator equivalent to $\\widehat{14}$')
(132, '', '2.44', 'Operator equivalent to $\\hat{7}$')
(133, '46', '2.58', 'Hadamard product $I \\ast \\iota$')
(134, '245', '2.67', 'This is operator "2.67" from ...')
(144, '~77,~78,~97', '2.43', 'Operator equivalent to $\\hat{6}$')
(149, '', '2.40', 'Operator equivalent to $\\hat{2}$ of AESZ.')
(176, '205', '2.69', 'This is operator "2.69" from ...')
(178, '182', '2.64', 'This is operator "2.64" from ...')
(184, '26', '2.61', 'A-incarnation: $X(1,1,1,1,2) \\subset  Grass(2,6)$')
(194, '140', '2.27', 'Hadamard product $D \\ast g$')
(210, '133', '2.20', 'Hadamard product A*f')
(217, '45', '2.1', 'Hadamard product $A \\ast a$, where $A$ is (:case 2.1.1)')
(233, '29', '2.53', 'Hadamard product $I \\ast \\gamma$')
(245, '58', '2.9', 'Hadamard product A*c')
(246, '36', '2.13', 'A*d')
(262, '134', '2.21', 'Hadamard product $B\\ast f$')
(266, '136', '2.23', 'Hadamard product $D \\ast f$')
(279, '41', '2.54', 'Hadamard product $I \\ast \\delta$')
(280, '184\n', '2.57', 'Hadamard product $I \\ast \\eta$')
(291, '137', '2.24', 'Hadamard product $A \\ast g$.')
(292, '', '2.70', 'Operator equivalent to (:aesz 255)')
(293, '25', '2.5', 'Hadamard product $A\\ast b$')
(301, '15', '2.2', 'Hadamard product $B\\ast a$.')
(320, '48', '2.14', 'B*d')
(322, '65', '2.16', 'Hadamard product D*d')
(343, '70', '2.10', 'Hadamard product $B\\ast c$.')
(344, '138', '2.25', 'Hadamard product $B\\ast g$.')
(371, '111', '2.17', 'B-Incarnations:')
(388, '135', '2.22', '>console.log("JSON Response:")')
(389, '24', '2.6', 'Hadamard product B*b')
(390, '110', '2.18', '>console.log("JSON Response:")')
(397, '62', '2.4', 'Hadamard product D*a')
(411, '18', '2.60', 'A-Incarnation: (1,1) and (2,2) intersection in $P^3 \\times P^3$')
(414, '16', '2.52', 'Hadamard product $I \\ast \\alpha$')
(422, '64', '2.12', 'Hadamard product $B\\ast c$.')
(423, '139', '2.26', 'Hadamard product $C \\ast g$')
(426, '142', '2.29', 'Hadamard product $A \\ast$')
(435, '49', '2.28', 'Hadamard product $A \\ast$')
(437, '38', '2.15', 'Hadamard product $C\\ast d$')
(438, '~67', '2.35', 'This is operator "2.35" from ...')
(446, '28', '2.62', 'A-incarnation: $X(1, 1, 1, 1, 1, 1) \\subset  Grass(3, 6)$')
(454, '68', '2.3', 'C*a')
(471, '185\n', '2.56', 'Hadamard product $I \\ast \\zeta$')
(477, '', '2.32', 'This is operator "2.32" from ...')
(484, '69', '2.11', 'Hadamard product $C\\ast c$')
(495, '63', '2.8', 'Hadamard product D*b')
(515, '42', '2.55', 'Hadamard product $I \\ast \\epsilon$')
(542, '51', '2.7', 'Hadamard product C*b')
(548, '7**', '2.31', 'This operator replaces AESZ 43, to which it is equivalent.')
'''

TryList = [7, 8, 17, 27, 38, 45, 54, 66, 68, 73, 75, 
           81, 82, 98, 100, 102, 104, 106, 107, 113, 
           126, 131, 132, 133, 134, 144, 149, 176, 178, 
           184, 194, 210, 217, 233, 245, 246, 262, 266, 
           279, 280, 291, 292, 293, 301, 320, 322, 343, 
           344, 371, 388, 389, 390, 397, 411, 414, 422, 
           423, 426, 435, 437, 438, 446, 454, 471, 477, 
           484, 495, 515, 542, 548]

lenTryList = 70

outlist1 = []
outlist2 = []
for n in range(lenTryList):
    AESZnum = TryList[n]
    thepfo = AESZ.AESZ(AESZnum).pfo
    thebound = min(abs(_) for _ in sympy.solve(thepfo.discriminant))
    solinfo = thepfo.all_sol(100)[0]
    possiblez = []
    for a in range(-30, 9):
        for b in range(-24, 9):
            for c in range(-6, 4):
                for d in range(-2, 3):
                    z0 = (sympy.Integer(2) ** a) * (sympy.Integer(3) ** b) * (sympy.Integer(5) ** c) * (sympy.Integer(7) ** d)
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

with open('solin.txt', 'w') as f:
    f.write("sollist = " + st1 + ";\n" + "zlist = " + st2 + ";")
