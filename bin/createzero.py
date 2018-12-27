import sys
import numpy as np

f = open(sys.argv[1])


x1 = []
x2 = []
x3 = []
x4 = []
x5 = []


for line in f:
    x1.append(line.split()[0])
    x2.append(line.split()[1])
    x3.append(line.split()[2])
    x4.append(line.split()[3])
    x5.append(line.split()[4])

length = len(x1)
y1 = np.zeros((length), dtype=np.int)

for i in range(length):
    y1[i]=1

y = []


f1 = open(sys.argv[2],"w")
for i in range(len(y1)):
    f1.write(str(y1[i])+"\n")
f1.close()