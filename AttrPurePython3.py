import AESZ
import PFTool
import sympy
import mpmath
from LLL import reduction

z2_1 = sympy.sqrt(2) - sympy.Integer(1)
z3_1 = sympy.Integer(2) - sympy.sqrt(3)
z5_1 = (sympy.sqrt(5) - sympy.Integer(1)) / 2
z7_1 = sympy.Integer(8) - 3 * sympy.sqrt(7)
z11_1 = sympy.Integer(10) - 3 * sympy.sqrt(11)
z13_1 = (sympy.sqrt(13) - sympy.Integer(3)) / 2
z17_1 = sympy.sqrt(17) - sympy.Integer(4)
z19_1 = sympy.Integer(170) - 39 * sympy.sqrt(19)

z2_2 = 1 / sympy.sqrt(2)
z3_2 = 1 / (sympy.sqrt(3) + sympy.Integer(1))
z5_2 = 1 / sympy.Integer(2)
z7_2 = 1 / (sympy.sqrt(7) + sympy.Integer(3))
z11_2 = 1 / (sympy.sqrt(11) + sympy.Integer(3))
z13_2 = 1 / sympy.Integer(2)
z17_2 = 2 / (sympy.sqrt(17) + sympy.Integer(3))
z19_2 = 1 / (sympy.sqrt(19) * 3 + sympy.Integer(13))

z2_3 = 1 / sympy.Integer(3)
z3_3 = 1 / sympy.sqrt(3)
z5_3 = 1 / sympy.Integer(3)
z7_3 = 1 / (sympy.sqrt(7) + sympy.Integer(2))
z11_3 = 1 / sympy.Integer(3)
z13_3 = 2 / (sympy.sqrt(13) + sympy.Integer(1))
z17_3 = 1 / sympy.Integer(3)
z19_3 = 1 / (sympy.sqrt(19) + sympy.Integer(4))

list1 = [z2_1, z3_1, z5_1, z7_1, z11_1, z13_1, z17_1, z19_1]
list2 = [z2_2, z3_2, z5_2, z7_2, z11_2, z13_2, z17_2, z19_2]
list3 = [z2_3, z3_3, z5_3, z7_3, z11_3, z13_3, z17_3, z19_3]
list4 = [['sqrt(2)'], ['sqrt(3)'], ['sqrt(5)'], ['sqrt(7)'], ['sqrt(11)'], ['sqrt(13)'], ['sqrt(17)'], ['sqrt(19)']]
list5 = [mpmath.sqrt(2), mpmath.sqrt(3), mpmath.sqrt(5), mpmath.sqrt(7), mpmath.sqrt(11), mpmath.sqrt(13), mpmath.sqrt(17), mpmath.sqrt(19)]

def recognize(inval, primesqrt):
    largenum = 10 ** 30
    s = [[1, 0, 0, int(largenum * inval)], [0, 1, 0, int(largenum * primesqrt)], [0, 0, 1, largenum]]
    return reduction(s)

def generate_testlist(bound, ithprime):
    testlist1 = [list1[ithprime] ** i for i in range(-4, 7)] + [-list1[ithprime] ** i for i in range(-4, 7)]
    if list2[ithprime] == 1/sympy.sqrt(2):
        testlist2 = [list2[ithprime] ** i for i in range(0, 6)]
    elif list2[ithprime] == 1/sympy.Integer(2):
        testlist2 = [list2[ithprime] ** i for i in range(0, 4)]
    else:
        testlist2 = [list2[ithprime] ** i * sympy.Integer(2) ** -j for j in range(0, 4) for i in range(-j, 4)]
    if list3[ithprime] == 1/sympy.sqrt(3):
        testlist3 = [list3[ithprime] ** i for i in range(0, 6)]
    elif list3[ithprime] == 1/sympy.Integer(3):
        testlist3 = [list3[ithprime] ** i for i in range(0, 4)]
    else:
        testlist3 = [list3[ithprime] ** i * sympy.Integer(3) ** -j for j in range(0, 4) for i in range(-j, 4)]
    testlist0 = [a * b * c for a in testlist1 for b in testlist2 for c in testlist3]
    testlist = [sympy.expand(x) for x in testlist0 if bound / 300 <= abs(PFTool.toc(x)) <= bound / 2]
    return testlist


for num in range(31, 552):
    pfo = AESZ.AESZ(num).pfo
    sol = pfo.mp_all_sol(230)
    bnd = abs(sol[0][0][-1] ** (-1/199))
    print('PFO = ', num, ', BND = ', bnd)
    with open('valid_results.txt', 'a') as f:
        f.write('PFO = ' + str(num) + ', BND ='+ str(bnd) + '\n')
    for i in range(8):
        testlist = generate_testlist(bnd, i)
        print(i, len(testlist))
        for z0 in testlist:
            val = pfo.attrval_MUM(200, z0)
            idt = recognize(val, list5[i])[0]
            if max(abs(x) for x in idt) < 500000:
                print('ATTR!', num, z0, val, idt)
                with open('valid_results.txt', 'a') as f:
                    f.write('ATTR! ' + str(num) + ',' + str(z0) + ',' + str(val) + ',' + str(idt) + '\n')

                    