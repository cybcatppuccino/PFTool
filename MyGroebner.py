"""Groebner bases algorithms. """


from __future__ import print_function, division

from sympy.polys.monomials import monomial_mul, monomial_lcm, monomial_divides, term_div
from sympy.polys.orderings import lex
from sympy.polys.polyerrors import DomainError

from sympy.polys.polytools import parallel_poly_from_expr
from sympy.polys.rings import PolyRing
from sympy.polys.polytools import Poly
from sympy import Integer

class PolyCombination:
    
    def __init__(self, inpoly=None, inlst=None):
        self.poly = inpoly
        self.lst = inlst
        
    def setsubdata(self):
        # print(self.lst)
        # self.rem = self.poly.rem
        self.LT = self.poly.LT
        self.LM = self.poly.LM
        self.iterterms = self.poly.iterterms
    
    def setdata(self, inpoly, num, all_num, ring):
        self.poly = ring.from_dict(inpoly.rep.to_dict())
        self.lst = [ring.zero for _ in range(all_num)]
        self.lst[num] = ring.one
        # print(self.poly, self.lst)
        # self.lst = [0 for _ in range(all_num)]
        # self.lst[num] = 1
        
        self.setsubdata()
        return self
    
    def rem_two(self, polycomb):
        outPoly = PolyCombination(self.poly.rem(polycomb.poly))
        mul = self.poly.quo(polycomb.poly)
        outPoly.lst = [self.lst[num] - polycomb.lst[num] * mul for num in range(len(self.lst))]
        # print("rem_two", self.poly, polycomb.poly)
        outPoly.setsubdata()
        return outPoly
    
    def rem(self, inlst):
        if len(inlst) == 0:
            return self
        elif len(inlst) == 1:
            return self.rem_two(inlst[0])
        else:
            return self.rem_two(inlst[0]).rem(inlst[1:])
    
    def monic(self):
        # return self.poly.monic()
        outPoly = PolyCombination(self.poly.monic())
        outPoly.lst = [term.quo_ground(self.poly.LC) for term in self.lst]
        # outPoly.lst = [term / self.poly.LC for term in self.lst]
        # print("monic")
        outPoly.setsubdata()
        return outPoly
    
    def mul_monom(self, inp):
        # return self.poly.mul_monom(inp)
        outPoly = PolyCombination(self.poly.mul_monom(inp))
        outPoly.lst = [term.mul_monom(inp) for term in self.lst]
        # print("mul_monom")
        outPoly.setsubdata()
        return outPoly
    
    # def notzero(self):
    #     return not (self.poly)

def MyMinus(polycomb1, polycomb2):
    outPoly = PolyCombination(polycomb1.poly - polycomb2.poly, [])
    outPoly.lst = [polycomb1.lst[num] - polycomb2.lst[num] for num in range(len(polycomb1.lst))]
    outPoly.setsubdata()
    return outPoly

def GroebnerBasis(F, *gens, **args):
    polys, opt = parallel_poly_from_expr(F, *gens, **args)
    ring = PolyRing(opt.gens, opt.domain, opt.order)
    polys = [PolyCombination().setdata(polys[num], num, len(polys), ring) for num in range(len(polys))]
    G = groebner(polys, ring, method='exbuchberger')
    Gaux = list(list(term.as_expr() for term in g.lst) for g in G)
    G = [g.poly.as_expr() for g in G]
    return G, Gaux


def groebner(seq, ring, method='exbuchberger'):

    _groebner_methods = {
        'exbuchberger': _exbuchberger,
    }

    _groebner = _groebner_methods[method]

    domain, orig = ring.domain, None

    if not domain.is_Field or not domain.has_assoc_Field:
        try:
            orig, ring = ring, ring.clone(domain=domain.get_field())
        except DomainError:
            raise DomainError("can't compute a Groebner basis over %s" % domain)
        else:
            seq = [ s.set_ring(ring) for s in seq ]

    G = _groebner(seq, ring)

    if orig is not None:
        G = [ g.clear_denoms()[1].set_ring(orig) for g in G ]

    return G

def _exbuchberger(f, ring):
    """
    Computes Groebner basis for a set of polynomials in `K[X]`. The original algorithm is the non-extended version. 
    """
    
    order = ring.order

    monomial_mul = ring.monomial_mul
    monomial_div = ring.monomial_div
    monomial_lcm = ring.monomial_lcm

    def select(P):
        # normal selection strategy
        # select the pair with minimum LCM(LM(f), LM(g))
        pr = min(P, key=lambda pair: order(monomial_lcm(f[pair[0]].LM, f[pair[1]].LM)))
        return pr

    def normal(g, J):
        # h = g.rem([ f[j] for j in J ])
        # print("normal = ", [ f[j] for j in J ])
        
        '''
        h = g
        for j in J:
            h = h.rem_two(f[j])
        '''
        
        h = g.rem([ f[j] for j in J ])
        # print("h1", g.poly, h.poly, [ f[j].poly for j in J ], bool(h.poly))

        if not h.poly:
        # if h.notzero():
            return None
        else:
            h = h.monic()

            if not h in I:
                I[h] = len(f)
                f.append(h)

            return h.LM, I[h]

    def update(G, B, ih):
        # update G using the set of critical pairs B and h
        # [BW] page 230
        h = f[ih]
        mh = h.LM

        # filter new pairs (h, g), g in G
        C = G.copy()
        D = set()

        while C:
            # select a pair (h, g) by popping an element from C
            ig = C.pop()
            g = f[ig]
            mg = g.LM
            LCMhg = monomial_lcm(mh, mg)

            def lcm_divides(ip):
                # LCM(LM(h), LM(p)) divides LCM(LM(h), LM(g))
                m = monomial_lcm(mh, f[ip].LM)
                return monomial_div(LCMhg, m)

            # HT(h) and HT(g) disjoint: mh*mg == LCMhg
            if monomial_mul(mh, mg) == LCMhg or (
                not any(lcm_divides(ipx) for ipx in C) and
                    not any(lcm_divides(pr[1]) for pr in D)):
                D.add((ih, ig))

        E = set()

        while D:
            # select h, g from D (h the same as above)
            ih, ig = D.pop()
            mg = f[ig].LM
            LCMhg = monomial_lcm(mh, mg)

            if not monomial_mul(mh, mg) == LCMhg:
                E.add((ih, ig))

        # filter old pairs
        B_new = set()

        while B:
            # select g1, g2 from B (-> CP)
            ig1, ig2 = B.pop()
            mg1 = f[ig1].LM
            mg2 = f[ig2].LM
            LCM12 = monomial_lcm(mg1, mg2)

            # if HT(h) does not divide lcm(HT(g1), HT(g2))
            if not monomial_div(LCM12, mh) or \
                monomial_lcm(mg1, mh) == LCM12 or \
                    monomial_lcm(mg2, mh) == LCM12:
                B_new.add((ig1, ig2))

        B_new |= E

        # filter polynomials
        G_new = set()

        while G:
            ig = G.pop()
            mg = f[ig].LM

            if not monomial_div(mg, mh):
                G_new.add(ig)

        G_new.add(ih)

        return G_new, B_new
        # end of update ################################

    if not f:
        return []

    # replace f with a reduced list of initial polynomials; see [BW] page 203
    f1 = f[:]

    while True:
        f = f1[:]
        f1 = []

        for i in range(len(f)):
            p = f[i]
            r = p.rem(f[:i])
            
            # print("r", r.poly, list(term.poly for term in f[:i]), bool(r.poly))
            if r.poly:
            # if not r.notzero():
                f1.append(r.monic())

        if len(f) == len(f1):
            break

    I = {}            # ip = I[p]; p = f[ip]
    F = set()         # set of indices of polynomials
    G = set()         # set of indices of intermediate would-be Groebner basis
    CP = set()        # set of pairs of indices of critical pairs

    for i, h in enumerate(f):
        I[h] = i
        F.add(i)

    #####################################
    # algorithm GROEBNERNEWS2 in [BW] page 232

    while F:
        # select p with minimum monomial according to the monomial ordering
        h = min([f[x] for x in F], key=lambda f: order(f.LM))
        ih = I[h]
        F.remove(ih)
        G, CP = update(G, CP, ih)

    # count the number of critical pairs which reduce to zero
    reductions_to_zero = 0

    while CP:
        ig1, ig2 = select(CP)
        CP.remove((ig1, ig2))

        h = spoly(f[ig1], f[ig2], ring)
        # ordering divisors is on average more efficient [Cox] page 111
        G1 = sorted(G, key=lambda g: order(f[g].LM))
        ht = normal(h, G1)

        if ht:
            G, CP = update(G, CP, ht[1])
        else:
            reductions_to_zero += 1

    ######################################
    # now G is a Groebner basis; reduce it
    Gr = set()

    for ig in G:
        ht = normal(f[ig], G - set([ig]))

        if ht:
            Gr.add(ht[1])

    Gr = [f[ig] for ig in Gr]

    # order according to the monomial ordering
    Gr = sorted(Gr, key=lambda f: order(f.LM), reverse=True)

    return Gr

def spoly(p1, p2, ring):
    """
    Compute LCM(LM(p1), LM(p2))/LM(p1)*p1 - LCM(LM(p1), LM(p2))/LM(p2)*p2
    This is the S-poly provided p1 and p2 are monic
    """
    LM1 = p1.LM
    LM2 = p2.LM
    LCM12 = ring.monomial_lcm(LM1, LM2)
    m1 = ring.monomial_div(LCM12, LM1)
    m2 = ring.monomial_div(LCM12, LM2)
    s1 = p1.mul_monom(m1)
    s2 = p2.mul_monom(m2)
    s = MyMinus(s1, s2)
    return s

def red_groebner(G, ring):
    """
    Compute reduced Groebner basis, from BeckerWeispfenning93, p. 216

    Selects a subset of generators, that already generate the ideal
    and computes a reduced Groebner basis for them.
    """
    def reduction(P):
        """
        The actual reduction algorithm.
        """
        Q = []
        for i, p in enumerate(P):
            h = p.rem(P[:i] + P[i + 1:])
            # print("h2", type(h), h, bool(h))
            
            if h.poly:
            # if not h.notzero():
                Q.append(h)

        return [p.monic() for p in Q]

    F = G
    H = []

    while F:
        f0 = F.pop()

        if not any(monomial_divides(f.LM, f0.LM) for f in F + H):
            H.append(f0)

    # Becker, Weispfenning, p. 217: H is Groebner basis of the ideal generated by G.
    return reduction(H)


if __name__ == "__main__":
    import sympy

    # sympy.Symbols
    x1 = sympy.Symbol('x1')
    x2 = sympy.Symbol('x2')
    x3 = sympy.Symbol('x3')
    y1 = sympy.Symbol('y1')
    y2 = sympy.Symbol('y2')
    y3 = sympy.Symbol('y3')
    z = sympy.Symbol('z')
    gb = GroebnerBasis([3*x1**2 - 2*x1*x3*z - 2*x1*x3 + x3**2*z, 2*x2*x3, -x1**2*z - x1**2 + 2*x1*x3*z + x2**2], \
                       [x1, x2, x3], domain='QQ(z)', order='grlex')
    print(sympy.GroebnerBasis([3*x1**2 - 2*x1*x3*z - 2*x1*x3 + x3**2*z, 2*x2*x3, -x1**2*z - x1**2 + 2*x1*x3*z + x2**2], \
                       [x1, x2, x3], domain='QQ(z)', order='grlex'))
