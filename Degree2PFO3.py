import AESZ
import PFTool
import sympy
import sys
sys.set_int_max_str_digits(100000)

TheList = [184,217,301,397,411,454]

def get_transed_pfo(innum):
    X = AESZ.AESZ(innum)
    op = X.pfo
    val1 = min(sympy.solve(op.discriminant))
    val2 = max(sympy.solve(op.discriminant))
    transval = sympy.simplify((1/val1 + 1/val2)/2)

    op1 = op.translation_inf()
    op2 = op1.translation(transval)
    op3 = op2.translation_inf()
    newtform = sympy.simplify(op3.tform * (transval*PFTool.z+1) ** 2)
    print(newtform)
    op4 = PFTool.PFO(newtform)
    return op4, transval, sympy.simplify(2/(1/val2 - 1/val1))

# print(get_transed_pfo(7)[0])


def get_transed_val(inval, theval):
    return sympy.simplify(1/(1/theval - inval))

outlist1 = []
outlist2 = []
for n in TheList:
    thepfo, inval, thebound = get_transed_pfo(n)
    solinfo = thepfo.all_sol(300)[0]
    outlist1.append(solinfo)
    possiblez = []
    for a in range(-12, 1):
        for b in range(-6, 1):
            for c in range(-2, 1):
                for d in range(-1, 1):
                    z0 = (sympy.Integer(2) ** a) * (sympy.Integer(3) ** b) * (sympy.Integer(5) ** c) * (sympy.Integer(7) ** d)
                    if abs(get_transed_val(inval, z0)) < thebound:
                        possiblez.append(get_transed_val(inval, z0))
                    z0 = -(sympy.Integer(2) ** a) * (sympy.Integer(3) ** b) * (sympy.Integer(5) ** c) * (sympy.Integer(7) ** d)
                    if abs(get_transed_val(inval, z0)) < thebound:
                        possiblez.append(get_transed_val(inval, z0))
    outlist2.append(possiblez)
    print(str(n) + " done! " + str(len(possiblez)))

st1 = str(outlist1)
st1 = st1.replace('[', '{')
st1 = st1.replace(']', '}')

st2 = str(outlist2)
st2 = st2.replace('[', '{')
st2 = st2.replace(']', '}')

with open('soldegree2pfo3.txt', 'w') as f:
    f.write("sollist = " + st1 + ";\n" + "zlist = " + st2 + ";")

