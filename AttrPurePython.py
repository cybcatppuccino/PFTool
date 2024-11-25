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

thelist = [2,3,22,29,34,44,51,86,90,105,128,129,159,177,182,243,260,299,325,378,392,398,420,439,441,478,523,529]

import AESZ
import sympy
import mpmath

for num in thelist:
    pfo = AESZ.AESZ(num).pfo
    disc = pfo.discriminant
    solveinv = [1/_ for _ in sympy.solve(disc)]
    solveinv.sort()
    lcm = sympy.lcm([abs(_) for _ in solveinv])
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
    allfactors = [-_ for _ in allfactorspos] + allfactorspos
    allfactors.sort()

    print(num, solveinv, lcm, allfactors)


pfo = AESZ.AESZ(447).pfo
z0 = sympy.Integer(1)/8
#print(pfo.eval_attr_Wronskian_MMA(z0, pr=True))

print(pfo.eval_attr_LLL(z0, pr=True))