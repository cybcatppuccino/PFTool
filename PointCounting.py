from numba import njit

primes = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997,1009,1013,1019,1021,1031,1033]

@njit
def count_points(p):
    k=-1
    count = 0
    for a in range(1, p):
        for b in range(a+1, p):
            for c in range(b+1, p):
                p1 = (b * c) % p
                p2 = b + c
                p3 = (k * a * p1) % p
                p4 = (1 + a + p2) % p
                p5 = a * (p1 + p2)
                p6 = (p1 + p5) % p
                for d in range(c+1, p):
                    if ((p4+d)*(p1*a+p6*d)-p3*d) % p == 0:
                        count += 24
    for a in range(1, p):
        b=a
        for c in range(1, p):
            for d in range(1, p):
                if (a != c) and (a != d) and (c != d) and ((1+a+b+c+d)*(b*c*(a+d) + a*d*(b+c+b*c))-k*a*b*c*d) % p == 0:
                    count += 6
    for a in range(1, p):
        b=a
        c=a
        for d in range(1, p):
            if (a != d) and ((1+a+b+c+d)*(b*c*(a+d) + a*d*(b+c+b*c))-k*a*b*c*d) % p == 0:
                    count += 4
    for a in range(1, p):
        b=a
        for c in range(a+1, p):
            d=c
            if (a != c) and ((1+a+b+c+d)*(b*c*(a+d) + a*d*(b+c+b*c))-k*a*b*c*d) % p == 0:
                    count += 6
    for a in range(1, p):
        b=a
        c=a
        d=a
        if ((1+a+b+c+d)*(b*c*(a+d) + a*d*(b+c+b*c))-k*a*b*c*d) % p == 0:
                    count += 1
    return count

lst = []
for p in primes:
    if 0 < p < 200:
        print(p)
        lst.append(count_points(p))

print(lst)