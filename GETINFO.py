f = open( 'hrefs.txt', 'r', encoding = 'utf-8' )
instr = f.read()
rs = []
n = 0
while ('href' in instr):
    num = instr.find('href')
    num = num + 8
    num2 = num
    while instr[num2] != '"':
        num2 += 1
    rs.append('https://cydb.mathematik.uni-mainz.de/' + instr[num:num2].replace('amp;', ''))
    instr = instr[num2:]
    n += 1