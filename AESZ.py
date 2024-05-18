import sympy
import PFTool

def get_raw_info(num):
    f = open( 'hrefs.txt', 'r', encoding = 'utf-8' )
    outstr = f.read()
    f.close

class AESZ:
    def __init__(self, num):
        self.num = num
        self.AESZnum = None
        self.newnum = None
        self.ss = (0, 0)
        
        