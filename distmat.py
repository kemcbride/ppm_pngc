#!/usr/bin/env python
import sys
homefile = open(sys.argv[1],'r')
output = open(sys.argv[2],'w')

x = 0
l = []
homefile.readline()
for line in homefile:
    l.append(line.strip('\n'))
    x+=1
    if x==2:
        for item in l:
            output.write(item+' ')
            x=0
            l = []
        output.write('\n')
    else:
        pass
homefile.close()
output.close()
