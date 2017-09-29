# Testing some values against the "known" working matlab code
from vpmfuncs import *

m = 0
p = 0
t = 12.0/100
c = 1

x, y = NACAAirfoil(m,p,t,c)

for k in x:
    print(k)
