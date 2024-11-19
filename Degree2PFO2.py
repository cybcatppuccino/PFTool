import AESZ
import PFTool
import sympy
import sys
sys.set_int_max_str_digits(100000)

TheList = [176, 178, 194, 245, 246, 291, 320, 322, 343, 344, 414, 422, 423, 437, 446, 484]

def get_transed_pfo(innum):
    X = AESZ.AESZ(innum)
    op = X.pfo
    val = min(sympy.solve(op.discriminant))
    transval = (1/(2*val))

    op1 = op.translation_inf()
    op2 = op1.translation(transval)
    op3 = op2.translation_inf()
    newtform = sympy.simplify(op3.tform * (transval*PFTool.z+1) ** 2)
    op4 = PFTool.PFO(newtform)
    return op4, val

# print(get_transed_pfo(7)[0])


def get_transed_val(inval, theval):
    return sympy.simplify(-1/(1/theval - 1/(2*inval)))

possiblez = []
for a in range(-5, 6):
    for b in range(-3, 4):
        for c in range(-2, 3):
            for d in range(-1, 2):
                z0 = (sympy.Integer(2) ** a) * (sympy.Integer(3) ** b) * (sympy.Integer(5) ** c) * (sympy.Integer(7) ** d)
                if 1 / 20 <= z0 <= 20:
                    possiblez.append(z0)

outlist1 = []
outlist2 = []
for n in TheList:
    thepfo, inval = get_transed_pfo(n)
    solinfo = thepfo.all_sol(300)[0]
    outlist1.append(solinfo)
    outlist2.append([get_transed_val(inval, -_*inval) for _ in possiblez])
    print(str(n) + " done! " + str(len(possiblez)))

st1 = str(outlist1)
st1 = st1.replace('[', '{')
st1 = st1.replace(']', '}')

st2 = str(outlist2)
st2 = st2.replace('[', '{')
st2 = st2.replace(']', '}')

with open('soldegree2pfo2.txt', 'w') as f:
    f.write("sollist = " + st1 + ";\n" + "zlist = " + st2 + ";")

