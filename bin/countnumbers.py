import sys
import numpy as np

f = open(sys.argv[1])
f1 = open(sys.argv[2])


x1 = []
x2 = []
x3 = []
x4 = []

z1 = []
z2 = []
z3 = []
z4 = []
z5 = []


for line in f:
    x1.append(line.split()[0])
    x2.append(line.split()[1])
    x3.append(line.split()[2])
    x4.append(line.split()[3])


for line in f1:
    z1.append(line.split()[0])
    z2.append(line.split()[1])
    z3.append(line.split()[2])
    z4.append(line.split()[3])
    z5.append(line.split()[4])

def notexit(i,delt):
    if len(delt) == 0:
        return True
    else:
        for j in range(len(delt)):
            if i == delt[j]:
                return False
        return True

def chachong(x1,x2,x3,x4,z1,z2,z3,z4,z5):
    dele = []
    length = len(x1)
    total = length + 1
    lengthTarget = len(z1)
    y1 = np.zeros((total), dtype=np.int)
    sum = 0
    for i in range(length):
        for j in range(lengthTarget):
            if (abs(int(x2[i])-int(z2[j]))<5 and abs(int(x3[i])-int(z3[j]))<5) and x1[i] == z1[j] and x4[i] == z4[j]:
                y1[i] = 1
                sum = sum + 1
    y1[length] = sum
    return y1
y1 = chachong(x1,x2,x3,x4,z1,z2,z3,z4,z5)
y = []
#import os
#os.mknod("C:/Users/tgh/Documents/WeChat Files/Tao-GuiHua/Files/output.bed")


f2 = open(sys.argv[3],"w")
for i in range(len(y1)):
    f2.write(str(y1[i])+"\n")
f2.close()