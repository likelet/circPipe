import sys

f = open(sys.argv[1])


x1 = []
x2 = []
x3 = []
x4 = []
x5 = []
x6 = []


for line in f:
    x1.append(line.split()[0])
    x2.append(line.split()[1])
    x3.append(line.split()[2])
    x4.append(line.split()[3])
    x5.append(line.split()[4])
    x6.append(line.split()[5])
    
def notexit(i,delt):
    if len(delt) == 0:
        return True
    else:
        for j in range(len(delt)):
            if i == delt[j]:
                return False
        return True

def chachong(x1,x2,x3,x4,x5,x6):
    y1 = []
    y2 = []
    y3 = []
    y4 = []
    y5 = []
    y6 = []
    dele = []
    length = len(x1)
    for i in range(length):
        if notexit(i,dele) == True:
            same = []
            same.append(i)
            for j in range(length-1-i):
                if (abs(int(x2[i])-int(x2[j+1+i]))<5 and abs(int(x3[i])-int(x3[j+1+i]))<5) and x1[i] == x1[j+1+i] and x6[i] == x6[j+1+i]:
                    print('---------')
                    print('i=%d,j=%d',i,j)
                    print(x2[i])
                    print(x2[j+1+i])
                    print(abs(int(x2[i])-int(x2[j+1+i])))
                    #print(abs(int(x3[i])-int(x3[j+1+i])>5))为什么是0
                    #same.append(j)
                    same.append(j+1+i)
                    dele.append(j+1+i)
            if len(same)>1:
                max_same = x5[same[0]]
                max_id = 0
                for m in range(len(same)):
                    if x5[same[m]]>max_same:
                         max_same = x5[same[m]]
                         max_id = m
                y1.append(x1[same[max_id]])
                y2.append(x2[same[max_id]])
                y3.append(x3[same[max_id]])
                y4.append(x4[same[max_id]])
                y6.append(x6[same[max_id]])
                sum = 0
                for k in range(len(same)):
                    sum = sum + int (x5[same[k]])
                y5.append(sum)
            else:
                y1.append(x1[i])
                y2.append(x2[i])
                y3.append(x3[i])
                y4.append(x4[i])
                y6.append(x6[i])
                y5.append(x5[i])
    return y1,y2,y3,y4,y5,y6
y1,y2,y3,y4,y5,y6 = chachong(x1,x2,x3,x4,x5,x6)
y = []
#import os
#os.mknod("C:/Users/tgh/Documents/WeChat Files/Tao-GuiHua/Files/output.bed")


f1 = open(sys.argv[2],"w")
for i in range(len(y1)):
    f1.write(str(y1[i])+"\t"+str(y2[i])+"\t"+str(y3[i])+"\t"+str(y6[i])+"\t"+str(y5[i])+"\n")
f1.close()