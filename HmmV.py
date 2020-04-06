import sys
import os


from pip._vendor.distlib.compat import raw_input
import numpy

#Open and Read Fasta File, reconstruction input data
with open(sys.argv[2], 'r') as f2:
    line = f2.readline()
    title = ''
    sequence = ''
    while line:
        line = line.rstrip('\n')
        if '>' in line:
            title = line
        else:
            sequence = sequence + line;
        line = f2.readline()


#Open and Read HMM File
with open(sys.argv[1], 'r') as f1:
    lines = f1.readlines()


#Split HMM File
hmm_list =[]
for line in lines:
    line = line.strip('\n')
    #print(line)
    hmm_line = line.split()
    hmm_list.append( hmm_line)
#print(hmm_list)

#HMM row 1, setting state, symbol, symbol name list
state=int(hmm_list[0][0])
symbol=int(hmm_list[0][1])
#print(state,symbol)
Name={}
for i in range(len(hmm_list[0][2])):
    Name[hmm_list[0][2][i]]=i
#print(Name)

#Generate State name list
Q=[]
for i in range(state):
    Q.append(chr(i + 65))
#print(Q)

#HMM row 2 generate and intialize the Initial prob table
I=[]
for i in hmm_list[1]:
    I.append(float(i))
#print(I)

#Genarate 2 table , Transition and Emission
A=numpy.zeros(shape=(state,state),dtype=float)
E=numpy.zeros(shape=(state,symbol),dtype=float)

#HMM other rows, A,E table initialization
for i in range(2,len(hmm_list)):
    count = 0
    for j in hmm_list[i]:
        if count>state-1:
            E[i-2][count-state]=float(j)
            count=count+1
        else:
            A[i-2][count]=float(j)
            count=count+1
#print(A)
#print(E)

# use log2 version
Alog=numpy.log2(A)
Elog=numpy.log2(E)
Ilog=numpy.log2(I)



#viterbi algorithm in log2 version
def viterbiLog(x):

    xL = map(Name.get, x)
    x=list(xL)
    nrow, ncol = 2, len(list(x))
    mat = numpy.zeros(shape=(nrow, ncol), dtype=float)  # generate prob table
    matTb = numpy.zeros(shape=(nrow, ncol), dtype=int)  # generate path table for backtrace

    # Fill in prob table first column
    for i in range(0, nrow):
        mat[i,0] = Elog[i, x[0]] + Ilog[i]

    # Fill the prob table and path table
    for j in range(1, ncol):
        for i in range(0, nrow):
            ep = Elog[i, x[j]]
            mx, mxi = mat[0, j - 1] + Alog[0, i] + ep, 0
            for i2 in range(1, nrow):
                pr = mat[i2, j - 1] + Alog[i2, i] + ep
                if pr > mx:
                    mx, mxi = pr, i2
            mat[i, j], matTb[i, j] = mx, mxi

    # Find final state with maximal prob
    omx, omxi = mat[0, ncol - 1], 0
    for i in range(1, nrow):
        if mat[i, ncol - 1] > omx:
            omx, omxi = mat[i, ncol - 1], i

    # Backtrace path table
    i, p = omxi, [omxi]
    flag=matTb[i,ncol-1]
    listPrint = [flag]         #list of state changing
    listB = [ncol]             #The position of state changing
    for j in range(ncol - 1, 0, -1):

        i = matTb[i, j]
        if i != flag:
            listB.append(j)
            listPrint.append(i)
            flag = i
        p.append(i)

    #print and write
    a = map(lambda x: Q[x], listPrint[::-1])   #Reverse and State number to state Name
    listPrintA = list(a)
    listB = listB[::-1]           #reverse
    width = len(str(listB[len(listB) - 1]) + "  " + str(listB[len(listB) - 1]) + listPrintA[i])  #print standard width
    countB = 0                  #The number of continuous B state
    total = 0                   #The  number of all continuous states
    startPositn = 1
    with open('result.txt', 'w', encoding='utf-8') as fw:
        for i in range(len(listB)):
            if listPrintA[i] == 'B':
                countB = countB + 1
            print("{:>{}}".format(str(startPositn) + " " + str(listB[i]) + " " + listPrintA[i], width))  #print, right align
            print("{:>{}}".format(str(startPositn) + " " + str(listB[i]) + " " + listPrintA[i], width), file=fw) #write in the txt file
            startPositn = listB[i] + 1
            total = total + 1

        p = ''.join(map(lambda x: Q[x], p[::-1]))  #reverse and state number to path , The whole path.
        TotalB=0;
        for i in p:
            if i=='B':
                TotalB=TotalB+1
        print('Total Number of B:'+str(TotalB), file=fw)   #write in file
        print('Total Number of B:',TotalB)  #print Total Number of B
        print('Continuous State B:'+str(countB), file=fw)
        print('Continuous State B:', countB)  #print Number of Continuous State B
        print('Total Continuous State:'+str(total), file=fw)
        print('Total Continuous State:', total)  #print Number of Total Continuous State
        return  p     #  path



x=sequence.upper()   #make all data become uppercase latter
viterbiLog(x)
