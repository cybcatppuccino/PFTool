s = ' '
inputlist = []
while s != '':
    s = input()
    inputlist.append(s)
outstr = '0'
for _ in inputlist[:-1]:
    d = _.find(']')
    lst = eval(_[:d+1])
    outstr += '+' + _[d+1:] + '*X^' + str(lst[0]) + '*Y^' + str(lst[1])
    if lst[0] != lst[1]:
        outstr += '+' + _[d+1:] + '*Y^' + str(lst[0]) + '*X^' + str(lst[1])
print(outstr)
    