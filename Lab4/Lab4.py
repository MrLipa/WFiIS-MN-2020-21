import matplotlib.pyplot as plt
import numpy as np
import copy
import random
import math
import time

def print_matrix(A,nazwa='A',rount=1):
    if type(A[0])==int or type(A[0])==float:
        print(f"__________Wektor {nazwa}[]__________")
        for i in range(N):
            print(round(A[i],rount))
        print()
        print()
    if type(A[0])==list:
        print(f"__________Wektor {nazwa}[]__________")
        for i in range(N):
            round_to_whole = [round(num, rount) for num in A[i]]
            print(f'[{round_to_whole}]')
        print()
        print()


def set_matrix(A):
    delta=2*L/N
    x=[0]*N
    for i in range(N):
        x[i]=-L+(i+1)*delta
        A[i][i]=pow(delta,-2)+pow(x[i],2)/2
    for i in range(1,N):
        A[i][i-1]=A[i-1][i]=-1/pow(delta,2)/2


def gershgorina_metod(A):
    sum=[0]*len(A)
    n=[0]*len(A)
    for i in range(len(A)):
        for j in range(len(A)):
            sum[i]+=abs(A[i][j])
        sum[i]-=A[i][i]
        n[i]=A[i][i]
    return -max(sum)-max(n),max(sum)+max(n)


def euclides_norm(A):
    sum=0
    for i in range(N):
        sum+=A[i]*A[i]
    return math.sqrt(sum)


def eigenvalues_bisection_method(A,n,start,end,iter):
    val=end
    for j in range(iter):
        count = 0
        w = [0] * (N + 1)
        w[0] = 1
        w[1] = A[0][0] - val
        for i in range(2, len(A) + 1):
            w[i] = (A[i - 1][i - 1] - val) * w[i - 1] - pow(A[i - 1][i - 2], 2) * w[i - 2]
        for i in range(N):
            if w[i] * w[i + 1] < 0:
                count += 1
        if count>=n:
            val=(val+start)/2
        else:
            val=val+val-start
            start=(val+start)/2
    return start


def eigenvector_bisection_method(A,val):
    x=[0]*len(A)
    x[0]=1
    x[1]=(val-A[0][0])/A[0][1]
    for i in range(2,len(A)):
        x[i]=((val-A[i-1][i-1])*x[i-1]-A[i-2][i-1]*x[i-2])/A[i-1][i]
    e=euclides_norm(x)
    for i in range(len(x)):
        x[i]=x[i]/e
    return x


L=-5
N=50
iter=50
colors=['k','y','b','c','r','g','m']

A=[[0 for i in range(N)] for i in range(N)]
set_matrix(A)
print(gershgorina_metod(A))

for j in range(1,6):
    start, end = gershgorina_metod(A)
    x = [0] * N
    y = [0] * N
    for i in range(N):
        x[i] = i + 1
    lam=eigenvalues_bisection_method(A, j, start, end, iter)
    print(f'Wartośc własna lambda {j} = {lam}'.format(j,lam))
    y = eigenvector_bisection_method(A, lam)

    plt.plot(x, y, label="(i,x(i))", color=colors[j-1])

    plt.grid()
    plt.title("Kolejne wartości komórek wektora własnego")
    #plt.xlabel('i')
    #plt.ylabel('x(i)')

    plt.legend(framealpha=1, frameon=True)

    plt.show()
for j in range(1,6):
    start, end = gershgorina_metod(A)
    x = [0] * N
    y = [0] * N
    for i in range(N):
        x[i] = i + 1
    lam=eigenvalues_bisection_method(A, j, start, end, iter)
    print(f'Wartośc własna lambda {j} = {lam}'.format(j,lam))
    y = eigenvector_bisection_method(A, lam)
    plt.scatter(x,y,s=25,color=colors[j-1])
    plt.plot(x, y, label='lambda={} = {:.2f}'.format(j,lam), color=colors[j-1])

    plt.grid()
    plt.title("Kolejne wartości komórek 5 wektorów własnych")
    #plt.xlabel('i')
    #plt.ylabel('x(i)')

    plt.legend(framealpha=1, frameon=True)

plt.show()
