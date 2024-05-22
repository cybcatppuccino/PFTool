from PFTool import *
from AESZ import *
from sympy import *

z = Symbol('z')
'''
p = (17 + 35 *z + 86 *z**2 + 19 *z**3 + 87 *z**4 + 92 *z**5 + 91 *z**6 + 33 *z**7 + 
 34 *z**8 + 40 *z**9 + 41 *z**10 + 27 *z**11 + 92 *z**12 + 14 *z**13 + 49 *z**14 + 
 80 *z**15 + 16 *z**16 + z**17 + 16 *z**18 + 50 *z**19 + 17 *z**20 + 25 *z**21 + 
 53 *z**22 + 80 *z**23 + 63 *z**24 + 55 *z**25 + 80 *z**26 + 58 *z**27 + 
 97 *z**28 + 77 *z**29 + 14 *z**30) * z
q = (84 * z + 58 *z + 14 *z**2 + 80 *z**3 + 74 *z**4 + 97 *z**5 + 83 *z**6 + 76 *z**7 + 
 56 *z**8 + 24 *z**9 + 30 *z**10 + 72 *z**11 + 41 *z**12 + 90 *z**13 + 54 *z**14 + 
 63 *z**15 + 21 *z**16 + 35 *z**17 + 64 *z**18 + 8 *z**19 + 44 *z**20 + 49 *z**21 + 
 77 *z**22 + 84 *z**23 + 45 *z**24 + 55 *z**25 + 81 *z**26 + 90 *z**27 + 
 35 *z**28 + 8 *z**29 + 15 *z**30)
'''
def eGCD(a, b):
    prevx, x = 1, 0
    prevy, y = 0, 1
    while b:
        rs = divmod(a, b)
        x, prevx = prevx - rs[0]*x, x
        y, prevy = prevy - rs[0]*y, y
        a, b = b, rs[1]
    if a >= 0:
        return a, (prevx, prevy)
    else:
        return -a, (-prevx, -prevy)

def integer_to_p_adic(inint, inprime):
    outtuple = (inint, 0)
    while True:
        rs = divmod(outtuple[0], inprime)
        if rs[1]:
            return outtuple
        else:
            outtuple = (rs[0], outtuple[1] + 1)

s = 2**15 * 3**22
print("Come")
for _ in range(10000):
    eGCD(216516513,411313182)
print("Done")
