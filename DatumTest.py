import AESZ

for n in range(1, 552):
    a = AESZ.AESZ(n)
    if a.pfo.yukawa(8) != a.yukawa:
        print("Warning! Number ", n)
    if n%10 == 0:
        print(n)

for n in range(1, 552):
    a = AESZ.AESZ(n)
    if a.pfo.qcoord(8)[1] != a.qcoord:
        print("Warning! Number ", n)
    if n%10 == 0:
        print(n)