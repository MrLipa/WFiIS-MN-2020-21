import matplotlib.pyplot as plt
import numpy as np
import copy
import random
import math
import time
from numpy.linalg import inv
import sympy as sym
import scipy.misc
from math import sin
from sympy import *

def fun1(x):
    return math.exp(-x**2)
def fun2(x):
    return math.cos(x)
def fun3(x):
    return math.sin(x)
def derivative1(fun,x,d=1e-6,n=1):
    return scipy.misc.derivative(fun,x,dx=d,n=n)
def derivative2(gfg_exp,x0,n=1):
    x, y = symbols('x y')
    dif=gfg_exp
    for i in range(n):
        dif = str(diff(dif, x))
    x=x0
    return eval(dif)
def calculate_c(fun,k):
    k=abs(k)
    if fun=="exp(-x**2)":
        if k % 2 == 0:
            return (-1) ** (int(k/2)) / math.factorial(int(k/2))
        else:
            return 0
    elif fun=="cos(x)":
        if k % 2 == 0:
            return (-1) ** (int(k/2)) / math.factorial(k)
        else:
            return 0
    elif fun=="sin(x)":
        if k % 2 == 1:
            return (-1) ** (int((k-1)/2)) / math.factorial(k)
        else:
            return 0

def factor_Pade(fun,N,M):
    A=np.zeros((M,M))
    y=np.zeros((M,1))
    b=np.zeros((M+1,1))
    a=np.zeros((N+1,1))
    for i in range(M):
        k=N+i+1
        y[i]=-calculate_c(fun,k)
    for i in range(M):
        for j in range(M):
            k=N-M+1+j+i
            A[i,j] = calculate_c(fun,k)
    x=np.linalg.inv(A)@y
    b[0]=1
    for i in range(M):
        b[i+1]=x[M-1-i]
    for i in range(N+1):
        for j in range(i+1):
            k=i-j
            c = calculate_c(fun,k)
            a[i]+=b[j]*c
    return a,b
def approximation_Pade(fun,N,M,x):
    a,b=factor_Pade(fun,N,M)
    y=[]
    for j in x:
        P = 0
        Q = 0
        for i in range(N+1):
            P+=a[i]*j**i
        for i in range(M+1):
            Q+=b[i]*j**i
        y.append(P/Q)
    return y
# a,b=-5,5
# tab=[(2, 2), (4, 4), (6, 6), (2, 4), (2, 6), (2, 8),(14,14)]
# tab1=["\u2082,\u2082","\u2084,\u2084","\u2086,\u2086","\u2082,\u2084","\u2082,\u2086","\u2082,\u2088","\u2081\u2084,\u2081\u2084"]
# for k in range(len(tab)):
#     N,M=tab[k]
#     x = list(np.linspace(a, b, num=2000))
#     y = [fun1(i) for i in x]
#     y1=approximation_Pade("exp(-x**2)",N,M,x)
#     plt.plot(x,y,color='black',linewidth=3,label=r'Teoria f(x)=$e^{-x^2}$')
#     plt.plot(x,y1,color='orange',label=f"Aproksymacja R{tab1[k]}")
#     plt.xlabel('x')
#     plt.ylabel('f(x)')
#     plt.title("Wykres dla N={} i M={}".format(N,M))
#     plt.legend(framealpha=1, frameon=True)
#     plt.show()
# a,b=-5,5
# tab=[(2, 2), (4, 4), (6, 6)]
# tab1=["\u2082,\u2082","\u2084,\u2084","\u2086,\u2086"]
# for k in range(len(tab)):
#     N,M=tab[k]
#     x = list(np.linspace(a, b, num=100))
#     y = [fun2(i) for i in x]
#     y1=approximation_Pade("cos(x)",N,M,x)
#     plt.plot(x,y,color='black',linewidth=3,label="Teoria f(x)=cos(x)")
#     plt.plot(x,y1,color='orange',label=f"Aproksymacja R{tab1[k]}")
#     plt.xlabel('x')
#     plt.ylabel('f(x)')
#     plt.title("Wykres dla N={} i M={}".format(N,M))
#     plt.legend(framealpha=1, frameon=True)
#     plt.show()
# a,b=-2*math.pi,2*math.pi
# tab=[(3, 3), (5, 5), (7, 7)]
# tab1=["\u2083,\u2083","\u2085,\u2085","\u2087,\u2087"]
# for k in range(len(tab)):
#     N,M=tab[k]
#     x = list(np.linspace(a, b, num=100))
#     y = [fun3(i) for i in x]
#     y1=approximation_Pade("sin(x)",N,M,x)
#     plt.plot(x,y,color='black',linewidth=3,label="Teoria f(x)=sin(x)")
#     plt.plot(x,y1,color='orange',label=f"Aproksymacja R{tab1[k]}")
#     plt.xlabel('x')
#     plt.ylabel('f(x)')
#     plt.title("Wykres dla N={} i M={}".format(N,M))
#     plt.legend(framealpha=1, frameon=True)
#     plt.show()

x=np.arange(-2*math.pi,2*math.pi,0.1)
y=np.sin(x)
fig,ax=plt.subplots()
ax.plot(x,y)
unit = 0.5
plt.xlabel('x')
plt.ylabel('f(x)')
x_tick = np.arange(-2, 2+unit, unit)
x_label = [r"$" + format(r, ".2g")+ r"\pi$" for r in x_tick]
ax.set_xticks(x_tick*np.pi)
ax.set_xticklabels(x_label, fontsize=9)
plt.plot(x,[0]*len(x),color='black')
plt.plot([0,0,0],[0,1.2,-1.2],color='black')
plt.show()

x=np.arange(-5,5,0.1)
y=np.cos(x)
fig,ax=plt.subplots()
ax.plot(x,y)
unit = 0.5
plt.xlabel('x')
plt.ylabel('f(x)')
x_tick = np.arange(-2, 2+unit, unit)
x_label = [r"$" + format(r, ".2g")+ r"\pi$" for r in x_tick]
ax.set_xticks(x_tick*np.pi)
ax.set_xticklabels(x_label, fontsize=9)
plt.plot(x,[0]*len(x),color='black')
plt.plot([0,0,0],[0,1.2,-1.2],color='black')
plt.show()

x = list(np.linspace(-5, 5, num=100))
y = [math.exp(-i**2) for i in x]
plt.plot(x,y)
plt.xlabel('x')
plt.ylabel('f(x)')
plt.plot(x,[0]*len(x),color='black')
plt.plot([0,0,0],[0,1.2,-1.2],color='black')
plt.show()