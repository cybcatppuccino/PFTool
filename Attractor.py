import PFTool
import AESZ
import sympy

X = AESZ.AESZ(175)
op = X.pfo
loc = op.calclocalexp()
at = -sympy.Integer(1) / 7
opat = op.translation(at)