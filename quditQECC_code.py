import numpy.linalg
import math
import random
import numpy
from matplotlib import pyplot as plt 

##sets the dimension throughout
##always is a prime
q=3


##converts integer n to a q-ary value with pts digits
def inttoarr(n,pts):
    vals=numpy.zeros(pts)
    for k in range(pts):
        vals[k]=int(n/math.pow(q,pts-1-k))
        n=n-vals[k]*math.pow(q,pts-1-k)
    return vals

##returns the full symp group
##takes a stabilizer "generator" set c and prints the full group
##this does not require the generators to be linearly independent
##another function removes the dependencies
def sympgroup2(c):
    sympgroup=[]
    powers=numpy.zeros(len(c))
    for i in range(int(math.pow(q,len(c)))):
        powers=inttoarr(i,len(c))
        out=numpy.zeros(len(c[0]))
        for j in range(len(c)):
            temp=numpy.zeros(len(c[0]))
            for k in range(len(c[0])):
                temp[k]=powers[j]*c[j][k]
            out=numpy.add(out,temp)
            out=out%q
        sympgroup.append(out.tolist())
    return sympgroup

##verifies that the elements of c satisfy the requirements for stabilizer code
def verify(c):
    #returning 0 means it's valid [all commute], anything else is invalid
    valid=0
    for i in range(len(c)):
        for j in range(len(c)):
            valid=valid+comm(c[i],c[j])
    return valid

##computes the commutator of two symp paulis
def comm(p1,p2):
    if(not(len(p1)==len(p2))):
        print("dim error")
        return
    comm=0
    n=int(len(p1)/2)
    for i in range(n):
        comm=comm+(p1[i]*p2[i+n]-p1[i+n]*p2[i])
    return comm%q

##computes the commutator of two symp paulis over the integers
def comminf(p1,p2):
    if(not(len(p1)==len(p2))):
        print("dim error")
        return
    comm=0
    n=int(len(p1)/2)
    for i in range(n):
        comm=comm+(p1[i]*p2[i+n]-p1[i+n]*p2[i])
    return comm

##this will go through the group and check if any sum is repeated
##this makes the "generators" c into a set of independent generators
def makegens2(c):
    powers=numpy.zeros(len(c))
    newc=c.copy()
    removed=numpy.zeros(len(c))
    for i in range(int(math.pow(q,len(c)))):
        powers=inttoarr(i,len(c))
        out=numpy.zeros(len(c[0]))
        tally=0
        for l in range(len(c)):
            tally=tally+((powers[l]>0) and removed[l])
        if(tally >0):
            continue
        for j in range(len(c)):
            temp=numpy.zeros(len(c[0]))
            for k in range(len(c[0])):
                temp[k]=powers[j]*c[j][k]
            out=numpy.add(out,temp)
        out=out%q
        out=out.astype(int)
        out=out.tolist()
        if((out in newc) and not(math.log(i,q)==int(math.log(i,q)))):
            newc.remove(out)
            removed[c.index(out)]=1
    return newc

##computes the pauli weight of a symp pauli
def pweight(s):
    weight=0
    n=int(len(s)/2)
    for k in range(n):
        weight=weight+((s[k]>0) or (s[k+n]>0))
    return weight

##returns the distance of the code through brute-force
##assumes gens are gens
##this method is slower, but will find the distance exactly no matter the size
def dist(c):
    dmin=len(c[0])
    symp=sympgroup2(c)
    for i in range(int(math.pow(q,len(c[0])))):
        if(math.log(i,q)==int(math.log(i,q))):
            print(i)
        err=inttoarr(i,len(c[0]))
        err=err.tolist()
        temp=0
        for j in range(len(c)):
            temp=temp+(comm(err,c[j])%q)
        ##if it has a nonzero syndrome, skip it
        if(not(temp==0)):
            continue
        tem=(err in symp)
        if(not(tem)):
            dmin=min(dmin,pweight(err))
    return dmin

##this returns the distance of the code by building up by Hamming weights
##it only goes up to d=3, for d>3 it just returns that the distance is higher
##for smaller distances this works much faster
def dist2(c):
    dmin=len(c[0])
    symp=sympgroup2(c)
    found=0
    weight=1
    while(not(found)):
        errs=hammstrs(int(len(c[0])/2),weight)
        if(errs==[]):
            found=1
            print("distance is larger than 3")
            return
        for i in range(len(errs)):
            temp=0
            for j in range(len(c)):
                temp=temp+(comm(errs[i],c[j])%q)
            if(not(temp==0)):
                continue
            tem=(errs[i] in symp)
            if(not(tem)):
                dmin=weight
                found=1
                return dmin
        weight=weight+1
    return dmin

##generates all q-ary Hamming strings of length n with weight w
def hammstrs(n,w):
    #we will set t=q later
    t=q
    allstrs=[]
    if(w==0):
        temp=[]
        for j in range(2*n):
            temp.append(0)
        return temp
    if(w==1):
        for j in range(n):
            for q1 in range(t):
                for q2 in range(t):
                    temp=[]
                    for k in range(j):
                        temp.append(0)
                    temp.append(q1)
                    for k in range(n-j-1):
                        temp.append(0)
                    for k in range(j):
                        temp.append(0)
                    temp.append(q2)
                    for k in range(n-j-1):
                        temp.append(0)
                    allstrs.append(temp)
        return allstrs
    if(w==2):
        for j1 in range(n-1):
            for j2 in range(j1,n):
                for q1 in range(t):
                    for q2 in range(t):
                        for q3 in range(t):
                            for q4 in range(t):
                                temp=[]
                                for k in range(j1):
                                    temp.append(0)
                                temp.append(q1)
                                for k in range(j2-j1-1):
                                    temp.append(0)
                                temp.append(q3)
                                for k in range(n-j2-1):
                                    temp.append(0)
                                for k in range(j1):
                                    temp.append(0)
                                temp.append(q2)
                                for k in range(j2-j1-1):
                                    temp.append(0)
                                temp.append(q4)
                                for k in range(n-j2-1):
                                    temp.append(0)
                                if(len(temp)==2*n):
                                    allstrs.append(temp)
        return allstrs
    if(w==3):
        for j1 in range(n-2):
            for j2 in range(j1,n-1):
                for j3 in range(j2,n):
                    for q1 in range(t):
                        for q2 in range(t):
                            for q3 in range(t):
                                for q4 in range(t):
                                    for q5 in range(t):
                                        for q6 in range(t):
                                            temp=[]
                                            for k in range(j1):
                                                temp.append(0)
                                            temp.append(q1)
                                            for k in range(j2-j1-1):
                                                temp.append(0)
                                            temp.append(q3)
                                            for k in range(j3-j2-1):
                                                temp.append(0)
                                            temp.append(q5)
                                            for k in range(n-j3-1):
                                                temp.append(0)
                                            for k in range(j1):
                                                temp.append(0)
                                            temp.append(q2)
                                            for k in range(j2-j1-1):
                                                temp.append(0)
                                            temp.append(q4)
                                            for k in range(j3-j2-1):
                                                temp.append(0)
                                            temp.append(q6)
                                            for k in range(n-j3-1):
                                                temp.append(0)
                                            if(len(temp)==2*n):
                                                allstrs.append(temp)
        return allstrs                                    
    return []

##returns the multiplicative inverse of a number mod q
def minv(x,q):
    if(x==0):
        return -1
    for i in range(1,q):
        if(int((i*x)%q)==1):
            return i
    return "error"

##this performs the effect of conjugation by dft on column col
def dft(c,col):
    temp=0
    for i in range(len(c)):
        temp=c[i][col]
        c[i][col]=c[i][col+int(len(c[0])/2)]
        c[i][col+int(len(c[0])/2)]=(-1*temp)%q
    return c

##finite field reduced row echelon form--used for putting the code into canonical form
def ffrref(c):
    outc=numpy.zeros([len(c),len(c[0])])
    for i in range(len(c)):
        x=c[i][i]
        inv=minv(x,q)
        if(inv==-1):
            temp= outc.astype(int)
            return temp.tolist()
        for j in range(len(c[0])):
            outc[i][j]=(inv*c[i][j])%q
        for k in range(len(c)):
            if(not(k==i)):
                for j in range(len(c[0])):
                    outc[k][j]=(c[k][j]-c[k][i]*c[i][j])%q
    outc=outc.astype(int)
    return outc.tolist()

##this completes turning the code into canonical form
def symptocon(c):
    n = len(c[0]) // 2
    for i in range(len(c)):
        if((c[i][i] == 0)):
            stop=0
            for j in range(i, len(c)):
                if(stop>0):
                    continue
                f=n+i
                if(not(c[j][i] == 0)):
                    temp = c[i]
                    c[i] = c[j]
                    c[j] = temp
                    stop=1
                #need to do the Fourier Transform
                elif(not(c[j][f] == 0)):
                    temp = c[i]
                    c[i] = c[j]
                    c[j] = temp
                    c=dft(c, i)
                    stop=1
    return c

##generates a random stabilizer code over n qudits with k generators
def genrandoms(n,k):
    currcode=[]
    toadd=[]
    for m in range(2*n):
        toadd.append(random.randint(0,q-1))
    currcode.append(toadd)
    counter=0
    #makes sure the code doesn't get stuck in a bad seed forever
    while(len(currcode)<k and counter<50000):
        toadd=[]
        counter=counter+1
        #if(counter%100==0):
        #    print(counter)
        for m in range(2*n):
            toadd.append(random.randint(0,q-1))
        include=1
        for i in range(len(currcode)):
            if(not(comm(currcode[i],toadd)%q==0)):
               include=0
        if(include):
            temp=currcode.copy()
            temp.append(toadd)
            if(len(makegens2(temp))>len(currcode)):
                currcode.append(toadd)
    return currcode

##this verifies that the code is valid over the integers
def verifyinf(c):
    #returning 0 means it's valid [all commute], anything else is invalid
    valid=0
    for i in range(len(c)):
        for j in range(len(c)):
            valid=valid+comminf(c[i],c[j])
    return valid

##this generates the L matrix used to transform a code into LDI form
def lmatrix(c):
    l=numpy.zeros([len(c),len(c)])
    for i in range(len(c)):
        for j in range(len(c)):
            if(i >= j):
                n = comminf(c[i], c[j])
                l[i][j] = n
    l=l.astype(int)
    l=l.tolist()
    return l

##this takes in the code (in canonical form) and adds the L matrix
##finish c+[0_{kxk} 0_{kx(n-k)} | L_{kxk} 0_{kx(n-k)}]
##c[i][j+n]=L[i][j] 0<i,j<k
def invariant(c, l):
    n = len(c[0]) // 2
    k = len(c)
    for i in range(k):
        for j in range(k):
            c[i][j+n] += l[i][j]
    return c

##transforms the code c into LDI form, prints the LDI code and its distance
def quicktest(c):
    c=symptocon(c)
    c=ffrref(c)
    c=invariant(c,lmatrix(c))
    print(c)
    print(dist2(c))
    return

##takes in the sympgroup2 and returns the list of Pauli weight corresponding to it
def sympweights(c):
    out=[]
    for i in range(len(c)):
        out.append(pweight(c[i]))
    return out

def second(e):
    return e[1]

def minweightgens(c):
    k=len(c)
    group=sympgroup2(c)
    weights=sympweights(group)
    paired= list(map(list,zip(group,weights)))
    paired.sort(key=second)
    mingens=[paired[1]]
    index=1
    while(len(mingens)<k):
        currgens=[]
        for t in range(len(mingens)):
            currgens.append(mingens[t][0])
        partialgroup=sympgroup2(currgens)
        if(paired[index][0] not in partialgroup):
            mingens.append(paired[index])
        index=index+1
    print(mingens)
    return mingens[-1][1]

def forplot(n,k,d,levels):
    global q
    q=levels
    c=genrandoms(n,k)
    if(len(c)<k):
        #print("didn't find enough gens")
        return
    dq=dist2(c)
    #print(dq)
    if(dq<d):
        #print("distance was too low")
        return
    print(c)
    return

def findB(c):
    B=0
    for i in range(len(c)):
        for j in range(len(c[0])):
            if(abs(c[i][j])>B):
                B=abs(c[i][j])
    return B

def pnaive(n,k,d,q):
    B=(2+(n-k)*(q-1))*(q-1)
    return math.pow(B,2*(d-1))*math.pow(2*(d-1),(d-1))

##a set of stabilizer elements is an array of 2n-length arrays
##this defines the 9-qubit shor code
c9=[
[1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1]
]

##5 particle code--for qudit reasons
c5=[
[1, 0, 0, 1, 0, 0, 1, 1, 0, 0],
[0, 1, 0, 0, 1, 0, 0, 1, 1, 0],
[1, 0, 1, 0, 0, 0, 0, 0, 1, 1],
[0, 1, 0, 1, 0, 1, 0, 0, 0, 1]
]

##5 level code
c11=[[1,0,0,0,0,3,0,2,2,1,1,0,0,0,0,0,0,0,0,0,0,0],
     [0,1,0,0,0,1,2,1,0,3,2,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,1,0,0,4,1,3,1,1,4,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,1,0,4,3,3,3,0,1,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,1,1,3,0,2,2,1,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,3,0,2,2,1,1],
     [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,2,1,0,3,2],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,4,1,3,1,1,4],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,4,3,3,3,0,1],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,3,0,2,2,1]
     ]

c11three=[
     [1,0,0,0,0,2,2,2,1,0,1,0,0,0,0,0,0,0,0,0,0,0],
     [0,1,0,0,0,1,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,1,0,0,2,1,2,0,1,2,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,1,0,0,1,2,2,2,1,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,1,1,2,2,2,1,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,2,2,2,1,0,1],
     [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,1,1,1],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,2,1,2,0,1,2],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,2,2,2,1],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,2,2,2,1,0]
     ]

#quadratic residue 1 code
cqr1=[[1,0,0,0,0,2,2,1,2,0,1,0,0,0,0,0,0,0,0,0,0,0],
      [0,1,0,0,0,1,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0],
      [0,0,1,0,0,2,1,2,0,2,1,0,0,0,0,0,0,0,0,0,0,0],
      [0,0,0,1,0,0,1,1,1,2,1,0,0,0,0,0,0,0,0,0,0,0],
      [0,0,0,0,1,1,2,1,2,1,0,0,0,0,0,0,0,0,0,0,0,0],
      [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,2,2,1,2,0,1],
      [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,1,1,1],
      [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,2,1,2,0,2,1],
      [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,2,1],
      [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,2,1,2,1,0]]

##Euclidean Hamming, I think
#ceh=[[1,0,0,2,0,2,0,0,0,0,0,0],
#     [0,1,0,2,2,1,0,0,0,0,0,0],
#     [0,0,1,1,1,0,0,0,0,0,0,0],
#     [0,0,0,0,0,0,1,0,0,2,0,2],
#     [0,0,0,0,0,0,0,1,0,2,2,1],
#     [0,0,0,0,0,0,0,0,1,1,1,0]]

print("printing min weight stuff")
q=3

#7
bettertest=[[1, 2, 1, 2, 0, 0, 0, 0, 0, 1, 2, 1, 2, 2],
[2, 1, 2, 0, 1, 1, 0, 1, 0, 1, 0, 2, 1, 2],
[2, 0, 1, 2, 1, 0, 1, 0, 2, 0, 2, 2, 0, 2],
[1, 1, 0, 1, 2, 2, 2, 0, 0, 0, 0, 0, 2, 0]]
print(minweightgens(bettertest))

#10
lastone=[[0, 0, 0, 0, 0, 0, 2, 1, 2, 2, 0, 2, 0, 2, 1, 1, 0, 0, 2, 2],
[2, 0, 0, 2, 1, 0, 2, 0, 1, 0, 2, 2, 2, 0, 1, 0, 2, 2, 0, 2],
[2, 2, 1, 0, 0, 2, 2, 2, 0, 1, 2, 1, 0, 1, 0, 2, 2, 0, 1, 1],
[2, 0, 0, 0, 0, 1, 0, 0, 0, 2, 1, 2, 2, 2, 2, 1, 1, 1, 1, 0],
[0, 2, 0, 0, 1, 1, 2, 0, 0, 2, 1, 1, 1, 0, 1, 2, 2, 1, 0, 1]]
print(minweightgens(lastone))
