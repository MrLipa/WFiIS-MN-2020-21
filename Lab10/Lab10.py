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
    return math.log(x**5+3*x**2+x+9)
def fun2(x):
    return x**6
def Powell(fun,a,b,x,n,eps=1e-6):
    xm=[0]*n
    F1=[0]*n
    F2=[0]*n
    for i in range(n):
        x.sort(reverse=True)
        F1[i]=(fun(x[1])-fun(x[0]))/(x[1]-x[0])
        F2[i]=( (fun(x[2])-fun(x[1]))/(x[2]-x[1])-(fun(x[1])-fun(x[0]))/(x[1]-x[0]) ) / (x[2]-x[0])
        xm[i]=(x[0]+x[1])/2 - F1[i]/(2*F2[i])
        m=[]
        for j in x:
            m.append(abs(j-xm[i]))
        x[m.index(max(m))]=xm[i]
    ### druga wersja
    # xm = []
    # F1 = []
    # F2 = []
    # i = 0
    # while(True):
    #     x.sort(reverse=True)
    #     F1.append((fun(x[1])-fun(x[0]))/(x[1]-x[0]))
    #     F2.append(( (fun(x[2])-fun(x[1]))/(x[2]-x[1])-(fun(x[1])-fun(x[0]))/(x[1]-x[0]) ) / (x[2]-x[0]))
    #     xm.append((x[0]+x[1])/2 - F1[-1]/(2*F2[-1]))
    #     m=[]
    #     for j in x:
    #         m.append(abs(j-xm[i]))
    #     if (abs(x[m.index(max(m))]-xm[-1])<eps):
    #         break
    #     x[m.index(max(m))]=xm[-1]
    #     i+=1
    # del xm[0]
    # del F1[0]
    # del F2[0]
    return xm,F1,F2
a=-1.5
b=1
h=0.01
n=8
eps=1e-3


x=[-0.5+i*h for i in range(3)]
xm,F1,F2=Powell(fun1,a,b,x,n)

x1=list(np.linspace(a, b, num=100))
y1=[fun1(i) for i in x1]
plt.scatter(xm[-1],fun1(xm[-1]),color="red",label="ekstremum lokalne")
plt.scatter(xm[0],fun1(xm[0]),color="blue",label="punkt startowy")
plt.plot(x1,y1,color="black",label="wykres funkcji")
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title(f'Wykres f(x)=ln(x\u2075+3x\u00B2+x+9)\n k={len(xm)} iteracji')
plt.legend(framealpha=1, frameon=True)
plt.show()

plt.plot([i for i in range(1,n+1)],[xm[i] for i in range(n)],'-ok')
plt.xlabel('k')
plt.ylabel('x\u2098(k)')
plt.title("Wykres zalezności kolejnych x\u2098 od numeru iteracji")
plt.show()

plt.plot([i for i in range(1,n+1)],[F1[i] for i in range(n)],'-ok')
plt.xlabel('k')
plt.ylabel('F[x\u2081,x\u2082](k)')
plt.title("Wykres zalezności kolejnych ilorazem różnicowym pierwszego rzędu\n od numeru iteracji")
plt.show()

plt.plot([i for i in range(1,n+1)],[F2[i] for i in range(n)],'-ok')
plt.xlabel('k')
plt.ylabel('F[x\u2081,x\u2082,x\u2083](k)')
plt.title("Wykres zalezności kolejnych ilorazem różnicowym drugiego rzędu\n od numeru iteracji")
plt.show()

n=9
x=[-0.9+i*h for i in range(3)]
xm,F1,F2=Powell(fun1,a,b,x,n,eps)

x1=list(np.linspace(a, b, num=100))
y1=[fun1(i) for i in x1]
plt.scatter(xm[-1],fun1(xm[-1]),color="red",label="ekstremum lokalne")
plt.scatter(xm[0],fun1(xm[0]),color="blue",label="punkt startowy")
plt.plot(x1,y1,color="black",label="wykres funkcji")
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title(f'Wykres f(x)=ln(x\u2075+3x\u00B2+x+9)\n k={len(xm)} iteracji')
plt.legend(framealpha=1, frameon=True)
plt.show()

plt.plot([i for i in range(1,n+1)],[xm[i] for i in range(n)],'-ok')
plt.xlabel('k')
plt.ylabel('x\u2098(k)')
plt.title("Wykres zalezności kolejnych x\u2098 od numeru iteracji")
plt.show()

plt.plot([i for i in range(1,n+1)],[F1[i] for i in range(n)],'-ok')
plt.xlabel('k')
plt.ylabel('F[x\u2081,x\u2082](k)')
plt.title("Wykres zalezności kolejnych ilorazem różnicowym pierwszego rzędu\n od numeru iteracji")
plt.show()

plt.plot([i for i in range(1,n+1)],[F2[i] for i in range(n)],'-ok')
plt.xlabel('k')
plt.ylabel('F[x\u2081,x\u2082,x\u2083](k)')
plt.title("Wykres zalezności kolejnych ilorazem różnicowym drugiego rzędu\n od numeru iteracji")
plt.show()
a=-2.5
b=2.5
n=100
x=[1.5+i*h for i in range(3)]
xm,F1,F2=Powell(fun2,a,b,x,n,eps)

x1=list(np.linspace(a, b, num=100))
y1=[fun2(i) for i in x1]
plt.scatter(xm[-1],fun2(xm[-1]),color="red",label="ekstremum lokalne")
plt.scatter(xm[0],fun2(xm[0]),color="blue",label="punkt startowy")
plt.plot(x1,y1,color="black",label="wykres funkcji")
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title(f'Wykres f(x)=ln(x\u2075+3x\u00B2+x+9)\n k={len(xm)} iteracji')
plt.legend(framealpha=1, frameon=True)
plt.show()

plt.plot([i for i in range(1,n+1)],[xm[i] for i in range(n)],'-ok')
plt.xlabel('k')
plt.ylabel('x\u2098(k)')
plt.title("Wykres zalezności kolejnych x\u2098 od numeru iteracji")
plt.show()

plt.plot([i for i in range(1,n+1)],[F1[i] for i in range(n)],'-ok')
plt.xlabel('k')
plt.ylabel('F[x\u2081,x\u2082](k)')
plt.title("Wykres zalezności kolejnych ilorazem różnicowym pierwszego rzędu\n od numeru iteracji")
plt.show()

plt.plot([i for i in range(1,n+1)],[F2[i] for i in range(n)],'-ok')
plt.xlabel('k')
plt.ylabel('F[x\u2081,x\u2082,x\u2083](k)')
plt.title("Wykres zalezności kolejnych ilorazem różnicowym drugiego rzędu\n od numeru iteracji")
plt.show()
