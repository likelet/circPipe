import sys

f = open(sys.argv[1])

createVar = locals()
for line in f:
    colnums = len(line.split())

for i in range(colnums):
    createVar['a'+str(i)] = []
    createVar['b'+str(i)] = []

onlyreads = colnums - 4
for i in range(onlyreads):
    createVar['sum'+str(i)] = 0


f = open(sys.argv[1])
for line in f:
    for i in range(colnums):
        createVar['a'+str(i)].append(line.split()[i])

def notexit(i,delt):
    if len(delt) == 0:
        return True
    else:
        for j in range(len(delt)):
            if i == delt[j]:
                return False
        return True

rownums = len(createVar['a0'])
dele = []
for i in range(rownums):
    if notexit(i,dele) == True:
        same = []
        same.append(i)
        for j in range(rownums-1-i):
            if abs(int(createVar['a1'][i])-int(createVar['a1'][j+1+i]))<5 and abs(int(createVar['a2'][i])-int(createVar['a2'][j+1+i]))<5 and createVar['a0'][i] == createVar['a0'][j+1+i] and createVar['a3'][i] == createVar['a3'][j+1+i]:
                print('---------')
                print('i=%d,j=%d',i,j)
                print(createVar['a1'][i])
                print(createVar['a1'][j+1+i])
                print(abs(int(createVar['a1'][i])-int(createVar['a1'][j+1+i])))
                print(abs(int(createVar['a2'][i])-int(createVar['a2'][j+1+i])))
                same.append(j+1+i)
                dele.append(j+1+i)
        if len(same)>1:
            max_same = createVar['a4'][same[0]]
            max_id = 0
            for m in range(len(same)):
                if createVar['a4'][same[m]]>max_same:
                    max_same = createVar['a4'][same[m]]
                    max_id = m
            createVar['b0'].append(createVar['a0'][same[max_id]])
            createVar['b1'].append(createVar['a1'][same[max_id]])
            createVar['b2'].append(createVar['a2'][same[max_id]])
            createVar['b3'].append(createVar['a3'][same[max_id]])
            for z in range(onlyreads):
                createVar['sum'+str(z)] = 0
            lengthsame = len(same)
            for k in range(lengthsame):
                for z in range(onlyreads):
                    createVar['sum'+str(z)] = int (createVar['sum'+str(z)]) + int (createVar['a'+str(z+4)][same[k]])
            for z in range(onlyreads):
                createVar['sum'+str(z)] = int (createVar['sum'+str(z)]) / lengthsame
            for z in range(onlyreads):
                createVar['b'+str(z+4)].append(createVar['sum'+str(z)])
        else:
            for x in range(colnums):
                createVar['b'+str(x)].append(createVar['a'+str(x)][i])

newrow = len(createVar['b0'])
f1 = open(sys.argv[2],"w")
for i in range(newrow):
    f1.write(str(createVar['b0'][i])+"_"+str(createVar['b1'][i])+"_"+str(createVar['b2'][i])+"_"+str(createVar['b3'][i]))
    for j in range(onlyreads):
        f1.write("\t"+str(createVar['b'+str(j+4)][i]))
    f1.write("\n")
f1.close()

f2 = open(sys.argv[3],"w")
for i in range(newrow):
    f2.write(str(createVar['b0'][i])+"\t"+str(createVar['b1'][i])+"\t"+str(createVar['b2'][i])+"\t"+"."+"\t"+"."+"\t"+str(createVar['b3'][i])+"\n")
f2.close()
