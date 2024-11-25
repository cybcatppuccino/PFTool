'''
(2, '', '4.31', 'This is operator "4.31" from ...') (64*z - 1)**2*(108*z - 1)**2
(3, '60', '4.26', 'Sporadic YY-Operator') (16*z - 1)**2*(108*z - 1)**2
(22, '~37', '4.13', 'YY-Operator equivalent to AESZ 37$=C \\ast \\alpha ~tilde c \\ast i$') (256*z - 1)**2*(1024*z - 1)**2
(29, '356\n', '4.32', 'Sporadic YY-Operator.') (108*z - 1)**2*(128*z - 1)**2
(34, '', '4.35', 'Operator equivalent to 3.34, equivalent to') -(16*z + 1)*(128*z - 1)**3
(44, '276', '4.55', 'Sporadic Operator.') (110592*z - 1)**2*(331776*z - 1)**2
(51, '255', '4.29', 'Sporadic YY-Operator.') (256*z + 1)**4
(86, '189', '4.27', 'Sporadic YY-Operator') (4*z - 1)**2*(256*z - 1)**2
(90, '~39', '4.1', 'YY-Operator equivalent to AESZ 39=$A \\ast \\alpha$.') (64*z - 1)**2*(256*z - 1)**2
(105, '353', '4.71', 'Sporadic Operator, reducible to 3.33, so not a primary operator.') (16*z - 1)*(64*z - 1)**3
(128, '', '4.7', 'YY-operator equivalent to (:AESZ 50), $\\tilde B \\ast \\alpha$') (108*z - 1)**2*(432*z - 1)**2
(129, '254', '4.49', 'Sporadic Operator.') (4096*z - 1)*(6912*z - 1)*(15360*z - 1)**2
(159, '~66', '4.19', 'YY-Operator equivalent to $AESZ 66 =$D \\ast \\alpha \\tilde c \\ast j$') (1728*z - 1)**2*(6912*z - 1)**2
(177, '', '4.77', 'Sporadic Operator.') (16*z - 1)**2*(24*z - 1)*(25*z - 1)
(182, '244', '4.28', 'Sporadic YY-Operator') (4*z - 1)**2*(108*z + 1)**2
(243, '239', '4.47', 'Sporadic Operator.') -(432*z + 1)*(1728*z + 1)**2*(3456*z - 1)
(260, '', '4.70', 'Sporadic Operator. There is a second') (20736*z - 1)*(65536*z - 1)*(221184*z + 1)**2
(299, '237', '4.46', 'This is operator "4.46" from ...') (64*z - 1)*(96*z + 1)**2*(864*z - 1)
(325, '278', '4.57', 'Sporadic Operator.') (162*z + 1)**2*(432*z - 1)*(729*z - 1)
(378, '241', '4.48', 'Sporadic Operator.') -(64*z + 1)*(384*z - 1)**2*(1728*z - 1)
(392, '288', '4.60', 'Sporadic Operator.') (1536*z + 1)**2*(4096*z - 1)*(6912*z - 1)
(398, '264', '4.53', 'Sporadic Operator.') (4096*z - 1)*(27648*z - 1)*(43008*z + 1)**2
(420, '277', '4.56', 'Sporadic Operator. There is a second MUM-point hiding at') -(1024*z + 1)*(4096*z - 1)*(6144*z + 1)**2
(439, '258', '4.52', 'Sporadic Operator.') (256*z - 1)*(512*z + 1)**2*(1024*z - 1)
(441, '350', '4.69', 'Sporadic Operator. There is a second MUM-point hiding at infinity, corresponding to Operator AESZ 351/4.70') (24*z + 1)**2*(81*z - 1)*(256*z - 1)
(478, '', '4.76', 'Sporadic Operator.') -(27*z + 1)*(2592*z + 1)**2*(13824*z - 1)
(523, '233', '4.45', 'Sporadic Operator.') (192*z - 1)**2*(432*z - 1)*(512*z - 1)
(529, '289', '4.61', 'Sporadic Operator.') (256*z - 1)*(5120*z + 1)**2*(16384*z - 1)
'''

thelist = [44,90,105,128,129,159,177,260,299,325,378,392,398,439,441,523,529]

import AESZ
import PFTool
import sympy
import mpmath
from LLL import reduction

for num in thelist:
    pfo = AESZ.AESZ(num).pfo
    disc = pfo.discriminant
    solveinv = [1/_ for _ in sympy.solve(disc)]
    lcm = sympy.lcm([abs(_) for _ in solveinv])
    solveinv.append(0)
    solveinv.sort()
    factorint = sympy.factorint(lcm)
    def factors(indic):
        # return all factors of an integer by its prime factorization
        if len(indic) == 0:
            return [1]
        indicnew = indic.copy()
        p = list(indicnew.keys())[-1]
        e = indicnew.pop(p)
        f = factors(indicnew)
        outlst = []
        for i in range(int(e + 1)):
            outlst += [p ** i * _ for _ in f]
        return outlst
    allfactorspos = factors(factorint)
    # allfactors = [-_ for _ in allfactorspos] + allfactorspos
    allfactors = allfactorspos
    allfactors.sort()

    def invinterval(x, y):
        if x is None:
            return (None if y == 0 else 1/sympy.Integer(y), 0)
        if y is None:
            return (0, None if x == 0 else 1/sympy.Integer(x))
        return (None if y == 0 else 1/sympy.Integer(y), None if x == 0 else 1/sympy.Integer(x))

    for n in range(-1, len(solveinv)):
        piece = []
        for _ in allfactors:
            if n == -1 and _ < solveinv[0]:
                piece.append((1/sympy.Integer(_), invinterval(None, solveinv[0])))
            elif n == len(solveinv)-1 and _ > solveinv[-1]:
                piece.append((1/sympy.Integer(_), invinterval(solveinv[-1], None)))
            elif _ > solveinv[n] and _ < solveinv[n+1]:
                piece.append((1/sympy.Integer(_), invinterval(solveinv[n], solveinv[n+1])))
        piece.sort(key=lambda x: x[0])

        def LLLtest(mat):
            ls = list((mat ** -1) * mpmath.matrix([[1], [0], [0], [0]]))
            tpj = 2 * mpmath.pi * mpmath.j
            vals = [mpmath.zeta(3), tpj ** 3, tpj ** 2, tpj ** 1]
            ls2 = [mpmath.re(ls[0]*vals[0]), mpmath.re(ls[0]*vals[1]), mpmath.re(ls[1]*vals[2]), mpmath.re(ls[2]*vals[3]), mpmath.re(ls[3]), \
                mpmath.im(ls[0]*vals[0]), mpmath.im(ls[0]*vals[1]), mpmath.im(ls[1]*vals[2]), mpmath.im(ls[2]*vals[3]), mpmath.im(ls[3])]
            l = len(ls2) // 2
            llllst = [[int(a==b) + int(ls2[a] * 10**(mpmath.mp.dps - 12)) * int(b==l) + int(ls2[a + l] * 10**(mpmath.mp.dps - 12)) * int(b==l+1) for b in range(l+2)] for a in range(l)]
            return ls2, reduction(llllst)
        
        def genzlist(pt1, pt2):
            zlist = []
            znow = pt1[0]
            while znow < pt2[0]:
                zlist.append(znow)
                if pt1[1][0] is None:
                    znow += PFTool.approx((pt1[1][1] - znow) / 5)
                elif pt1[1][1] is None:
                    znow += PFTool.approx((znow - pt1[1][0]) / 4)
                else:
                    znow += PFTool.approx(min((pt1[1][1] - znow) / 5, (znow - pt1[1][0]) / 4))
            zlist.append(pt2[0])
            print(zlist)
            return zlist
        
        def is_attr(rst2):
            s1 = sum(abs(_) for _ in rst2[0]) < 10**5
            s5 = min(abs(_) for _ in rst2[0]) < 10**3
            return s1 and s5

        if len(piece) > 0:
            print(num, disc, piece)

            #'''

            mat = pfo.eval_Wronskian0(piece[0][0], pr=True)
            rst1, rst2 = LLLtest(mat)
            outstr = ""
            if is_attr(rst2):
                print("Attr!", num, piece[0][0], rst2)
                outstr = "Attr_"
            with open('results/' + outstr + str(num) + '_' + str(piece[0][0]).replace('/', 'd') + '.txt', 'w') as f:
                f.write(PFTool.mptomma(mat ** -1))
            
            for i in range(1, len(piece)):
                zlist = genzlist(piece[i-1], piece[i])
                mat = PFTool.hololist(pfo, zlist, 128) * mat
                rst1, rst2 = LLLtest(mat)
                outstr = ""
                if is_attr(rst2):
                    print("Attr!", num, piece[i][0], rst2)
                    outstr = "Attr_"
                with open('results/' + outstr + str(num) + '_' + str(piece[i][0]).replace('/', 'd') + 'r.txt', 'w') as f:
                    f.write(PFTool.mptomma(mat ** -1))

            #'''

'''

pfo = AESZ.AESZ(447).pfo
z0 = sympy.Integer(1)/8
#print(pfo.eval_attr_Wronskian_MMA(z0, pr=True))

print(pfo.eval_attr_LLL(z0, pr=True))

'''