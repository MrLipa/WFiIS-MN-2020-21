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


def print_equation_matrix(A,x,b,rount=1):
    print(f"__________A*x=b__________")
    print()
    temp=int(N/2)
    for i in range(N):
        round_to_whole = [round(num, rount) for num in A[i]]
        if i==temp:
            print(round_to_whole, '* [', round(x[i],rount), ']', " = ", '[', round(b[i],rount), ']')
        else:
            print(round_to_whole, '  [', round(x[i],rount), ']', "   ", '[', round(b[i],rount), ']')
    print()
    print()



def set_matrix(A,b,m):
    for i in range(N):
        for j in range(N):
            if abs(i-j)<=m:
                A[i][j]=1/(1+abs(i-j))
            else:
                A[i][j]=0
        b[i]=i+1


def euclides_norm(A):
    sum=0
    for i in range(N):
        sum+=A[i]*A[i]
    return math.sqrt(sum)


def calulate_steepest_descent(A,b,x,eps=1e-6):
    r = [0 for i in range(N)]
    rk = []
    xk = []
    while True:
        licz=0
        mian=0

        for i in range(N):
            temp=0
            for j in range(N):
                temp += A[i][j] * x[j]
            r[i] = b[i] - temp

        for i in range(N):
            temp=0
            for j in range(N):
                temp += A[i][j] * r[j]
            mian += temp * r[i]

        for i in range(N):
            licz += r[i] * r[i]

        a = licz / mian

        for i in range(N):
            x[i] = x[i] + a * r[i]

        rk.append(euclides_norm(r))
        xk.append(euclides_norm(x))

        if euclides_norm(r)<eps:
            break

    k = len(rk)
    return rk, xk, k


def calculate_conjugate_gradient(A,b,x,eps=pow(10,-3)):
    r = [0 for i in range(N)]
    v = [0 for i in range(N)]

    rk = []
    xk = []

    for i in range(N):
        temp=0
        for j in range(N):
            temp += A[i][j] * x[j]
        r[i] = b[i] - temp
        v[i] = b[i] - temp

    while euclides_norm(r) > eps:
        licz = 0
        mian = 0

        for i in range(N):
            temp=0
            for j in range(N):
                temp += A[i][j] * v[j]
            mian += temp * v[i]


        for i in range(N):
            licz += r[i] * r[i]

        a = licz / mian

        for i in range(N):
            x[i] = x[i] + a * v[i]

        licz = 0
        mian = 0
        r1=copy.deepcopy(r)

        for i in range(N):
            temp=0
            for j in range(N):
                temp += A[i][j] * v[j]
            r[i] = r[i] - a * temp

        for i in range(N):
            licz += r[i] * r[i]

        for i in range(N):
            mian += r1[i] * r1[i]

        b = licz / mian

        for i in range(N):
            v[i]=r[i]+b*v[i]

        rk.append(euclides_norm(r))
        xk.append(euclides_norm(x))
    k=len(rk)
    return rk, xk, k

def calculate_GJ(A, b, x):
    A1 = copy.deepcopy(A)
    b1 = copy.deepcopy(b)
    for i in range(len(A)):
        w = A1[i][i]
        b1[i] = b1[i] / w
        for j in range(len(A[0])):
            A1[i][j] = A1[i][j] / w
        for k in range(len(A)):
            w2 = A1[k][i]
            if k != i:
                b1[k] = b1[k] - b1[i] * w2
                for j in range(len(A[0])):
                    A1[k][j] = A1[k][j] - w2 * A1[i][j]
    for i in range(len(A)):
        x[i] = b1[i]

def set_matrix(A,b,m):
    for i in range(N):
        for j in range(N):
            if abs(i-j)<=m:
                A[i][j]=1/(1+abs(i-j))
            else:
                A[i][j]=0
        b[i]=i+1


N=100
eps=1e-3
A=[[0 for i in range(N)] for i in range(N)]
b=[0 for i in range(N)]

x=[0 for i in range(N)]

set_matrix(A,b,5)

start_time = time.time()
rk,xk,k=calulate_steepest_descent(A,b,x,eps)
print("Czas działania programu dla x=0: %s sekund" % (time.time() - start_time))
print("Liczba iteracji k:", k)


f=[i for i in range(k)]
fig,ax=plt.subplots()
ax.scatter(f,rk,label="(x,y)",s=10,color='black')
ax.plot(f,rk,label="(x,y)",linewidth=1,color=[0,0,0,1])
ax.set_yscale("log")
ax.set_xlabel("k (iteracja)")
ax.set_ylabel("||rk||2",color="black",fontsize=14)
ax2=ax.twinx()
ax2.scatter(f,xk,label="(x,y)",s=10,color='blue')
ax2.plot(f,xk,label="(x,y)",linewidth=1,color='blue')
ax2.set_ylabel("||xk||2",color="blue",fontsize=14)
plt.show()

#
# x=[0 for i in range(N)]
# eps=1e-6
# start_time = time.time()
# rk,xk,k=calulate_steepest_descent(A,b,x,eps)
# print("Czas działania programu: %s sekund" % (time.time() - start_time))
# print("Liczba iteracji k:", k)
#
# f=[i for i in range(k)]
# fig,ax=plt.subplots()
# ax.scatter(f,rk,label="(x,y)",s=10,color='black')
# ax.plot(f,rk,label="(x,y)",linewidth=1,color='black')
# ax.set_yscale("log")
# ax.set_xlabel("k (iteracja)")
# ax.set_ylabel("||rk||2",color="black",fontsize=14)
# ax2=ax.twinx()
# ax2.scatter(f,xk,label="(x,y)",s=10,color='blue')
# ax2.plot(f,xk,label="(x,y)",linewidth=1,color='blue')
# ax2.set_ylabel("||xk||2",color="blue",fontsize=14)
# plt.show()
#
#
# x=[1 for i in range(N)]
# eps=1e-6
# start_time = time.time()
# rk,xk,k=calulate_steepest_descent(A,b,x,eps)
# print("Czas działania programu: %s sekund" % (time.time() - start_time))
# print("Liczba iteracji k:", k)
#
# f=[i for i in range(k)]
# fig,ax=plt.subplots()
# ax.scatter(f,rk,label="(x,y)",s=10,color='black')
# ax.plot(f,rk,label="(x,y)",linewidth=1,color='black')
# ax.set_yscale("log")
# ax.set_xlabel("k (iteracja)")
# ax.set_ylabel("||rk||2",color="black",fontsize=14)
# ax2=ax.twinx()
# ax2.scatter(f,xk,label="(x,y)",s=10,color='blue')
# ax2.plot(f,xk,label="(x,y)",linewidth=1,color='blue')
# ax2.set_ylabel("||xk||2",color="blue",fontsize=14)
# plt.show()
#
#
# x=[0 for i in range(N)]
# eps=1e-6
# start_time = time.time()
# rk,xk,k=calculate_conjugate_gradient(A,b,x,eps)
# print("Czas działania programu: %s sekund" % (time.time() - start_time))
# print("Liczba iteracji k:", k)
#
# f=[i for i in range(k)]
# fig,ax=plt.subplots()
# ax.scatter(f,rk,label="(x,y)",s=10,color='black')
# ax.plot(f,rk,label="(x,y)",linewidth=1,color='black')
# ax.set_yscale("log")
# ax.set_xlabel("k (iteracja)")
# ax.set_ylabel("||rk||2",color="black",fontsize=14)
# ax2=ax.twinx()
# ax2.scatter(f,xk,label="(x,y)",s=10,color='blue')
# ax2.plot(f,xk,label="(x,y)",linewidth=1,color='blue')
# ax2.set_ylabel("||xk||2",color="blue",fontsize=14)
# plt.show()
#
#
# x=[0 for i in range(N)]
#
# start_time = time.time()
# rk,xk,k=calculate_GJ(A,b,x)
# print("Czas działania programu: %s sekund" % (time.time() - start_time))
