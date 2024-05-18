# To GetInfo From the WebPage:
# https://cydb.mathematik.uni-mainz.de/browse.php?superseeker

# https://cydb.mathematik.uni-mainz.de/?search=true&amp;m=lookup&amp;format=json&amp;superseeker=-1/3%2C-5/3&amp;searchButton=search
# https://cydb.mathematik.uni-mainz.de/?search=true&m=lookup&format=json&superseeker=-1/3%2C-5/3&searchButton=search

import requests
import GETINFO

rs = GETINFO.rs
print(rs[0])

for _ in range(len(rs)):
    r = requests.get(rs[_])
    file_name = '.\\info\\' + str(_ + 1) + '.txt'
    f = open(file_name, 'w')
    f.write(r.text)
    f.close()
    print(_ + 1, " finished!")