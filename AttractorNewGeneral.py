import AESZ
import sys
import PFTool
import sympy

sys.set_int_max_str_digits(100000)

num = 395
deg = 3

X = AESZ.AESZ(num)
op = X.pfo
print(op)
print(op.calclocalexp())
val = min(sympy.solve(op.discriminant))
transval = (1/(2*val))
print(val)

op1 = op.translation_inf()
op2 = op1.translation(transval)
op3 = op2.translation_inf()
newtform = sympy.simplify(op3.tform * (transval*PFTool.z+1) ** deg)
print(newtform)
op4 = PFTool.PFO(newtform)
print(op4)

solinfo = str(op4.all_sol(1200)[0])
solinfo = solinfo.replace('[', '{')
solinfo = solinfo.replace(']', '}')

with open('solingeneral.txt', 'w') as f:
    f.write('sol=' + str(solinfo) + ';')