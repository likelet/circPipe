import sys
import numpy as np

f = open(sys.argv[1])
f1 = open(sys.argv[2])


x1 = []

z1 = []
z2 = []


for line in f:
    x1.append(line.split()[0])


for line in f1:
    z1.append(line.split()[0])
    z2.append(line.split()[1])

def notexit(i,delt):
    if len(delt) == 0:
        return True
    else:
        for j in range(len(delt)):
            if i == delt[j]:
                return False
        return True

def chachong(x1,z1,z2):
    dele = []
    length = len(x1)
    lengthTarget = len(z1)
    y1 = np.zeros((length), dtype=np.int)
    for i in range(length):
        for j in range(lengthTarget):
            if x1[i] == z1[j]:
                y1[i] = z2[j]
    return y1
y1 = chachong(x1,z1,z2)
y = []
#import os
#os.mknod("C:/Users/tgh/Documents/WeChat Files/Tao-GuiHua/Files/output.bed")


f2 = open(sys.argv[3],"w")
for i in range(len(y1)):
    f2.write(str(y1[i])+"\n")
f2.close()