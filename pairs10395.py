N = 12 # also try N = 6
def generate_pairings(elements):
    if not elements:
        return [[]]
    first, pairs = elements[0], []
    for i in range(1, len(elements)):
        pair = (first, elements[i])
        remaining = elements[1:i] + elements[i+1:]
        for rest in generate_pairings(remaining):
            pairs.append([pair] + rest)
    return pairs

all_pairs = generate_pairings(list(range(N)))

def act2(pairs, num):
    for pair in pairs:
        if num in pair:
            return pair[0] + pair[1] - num
act3 = lambda num: 3*(num // 3) + ((num + 1) % 3)
act = lambda pairs, num: act3(act2(pairs, num))

def componentsize(pairs):
    s0, s1 = set([0]), None
    while s1 != s0:
        if s1 != None:
            s0 = s1
        s1 = s0.union(set(act3(x) for x in s0)).union(set(act2(pairs, x) for x in s0))
    return len(s0)

def cycletype(pairs):
    outlst = list()
    s = set(range(N))
    while len(s) > 0:
        p = s.pop()
        l, q = 1, act(pairs, p)
        while q in s:
            s.remove(q)
            l, q = l + 1, act(pairs, q)
        outlst.append(l)
        outlst.sort(key=lambda x:-x)
    return str(outlst)

countdict = dict()
for pairs in all_pairs:
    if componentsize(pairs) == N:
        c = cycletype(pairs)
        countdict[c] = 1 if c not in countdict else countdict[c] + 1
print(f"Total pairings: {len(all_pairs)}")
print(countdict)

'''
Total pairings: 10395
{'[5, 5, 1, 1]': 972, '[8, 2, 1, 1]': 972, '[9, 1, 1, 1]': 648,
'[11, 1]': 1944, '[6, 3, 2, 1]': 1944, '[6, 6]': 486, '[10, 2]': 972,
'[4, 4, 2, 2]': 486, '[8, 4]': 486, '[9, 3]': 648, '[3, 3, 3, 3]': 162}
'''


