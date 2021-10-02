import matplotlib.pyplot as plt
import numpy as np
import copy
import random
import math
import time
from numpy.linalg import inv
import sympy as sym

def fun1(x):
    return 1/(1+x**2)

def fun2(x):
    return math.cos(2*x)

def derivative_second(x,fun,dx=0.01):
    return (fun(x-dx)-2*fun(x)+fun(x+dx))/dx**2

def wyzM(knots,fun):
    h=[]
    for i in range(1,len(knots)):
        h.append(knots[i]-knots[i-1])
    lam=[]
    ni=[]
    for i in range(len(h)-1):
        lam.append(h[i+1]/(h[i]+h[i+1]))
        ni.append(1-h[i+1]/(h[i]+h[i+1]))
    d=[0]*len(knots)
    for i in range(1,len(knots)-1):
        d[i]=( 6/(h[i-1]+h[i])*((fun(knots[i+1])-fun(knots[i]))/h[i] - (fun(knots[i])-fun(knots[i-1]))/h[i-1] ))

    A=np.eye(len(knots))
    for i in range(len(knots)):
        A[i,i]+=1
    for i in range(len(knots)-2):
        A[i+1,i+2]=lam[i]
    for i in range(len(knots)-2):
        A[i+1,i]=ni[i]

    d1=[0]*len(knots)
    for i in range(1,len(knots)-1):
        d1[i]=3
    m=np.zeros((1,len(knots)))
    m=d@np.linalg.inv(A)
    return m
def wyzSx(knots,fun,x):
    i=0
    m=wyzM(knots,fun)
    h = []
    for i in range(1, len(knots)):
        h.append(knots[i] - knots[i - 1])

    for k in range(len(knots)-1):
        if knots[k]<x<knots[k+1]:
            i=k+1
        if x<=min(knots):
            i=1
        if x>=max(knots):
            i=len(knots)-1
    A=(fun(knots[i])-fun(knots[i-1]))/h[i-1]-h[i-1]/6*(m[i]-m[i-1])
    B=fun(knots[i-1])-m[i-1]*h[i-1]**2/6
    return m[i-1]*(knots[i]-x)**3/(6*h[i-1]) + m[i]*(x-knots[i-1])**3/(6*h[i-1]) + A*(x-knots[i-1])+B
def mistake(y1,y2):
    e=0
    for i in range(len(y1)):
        if(abs(y1[i]-y2[i])>e):
            e=abs(y1[i]-y2[i])
    return e

n = [5, 8, 15,21]
a = -5
b = 5
i=1
for j in n:
    knots = np.linspace(a, b, num=j)
    x = list(np.linspace(a, b, num=100))
    y1 = [wyzSx(knots,fun1, i) for i in x]
    plt.plot(x, y1, color='blue',label="Interpolacja")
    plt.scatter(knots, [fun1(i) for i in knots], color='red', s=40,label="węzły interpolacyjne")
    y2 = [fun1(x[i]) for i in range(len(x))]
    plt.plot(x, y2, color='black',label=r'f(x)=$\frac{1}{1+x^2}$')
    plt.grid()
    plt.title(f"Wykres {i}: Węzły równoodległe n={j},       h={round(knots[1]-knots[0],2)}      \u03B5\u2098\u2090\u2093={round(mistake(y1,y2),2)}")
    i+=1
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.legend(framealpha=1, frameon=True)
    plt.show()
for j in n:
    knots = np.linspace(a, b, num=j)
    x = list(np.linspace(a, b, num=100))
    y1 = [wyzSx(knots,fun2, i) for i in x]
    plt.plot(x, y1, color='blue',label="Interpolacja")
    plt.scatter(knots, [fun2(i) for i in knots], color='red', s=40,label="węzły interpolacyjne")
    y2 = [fun2(x[i]) for i in range(len(x))]
    plt.plot(x, y2, color='black',label=r'f(x)=cos(2x)$')
    plt.grid()
    plt.title(f"Wykres {i}: Węzły równoodległe n={j},       h={round(knots[1] - knots[0], 2)}      \u03B5\u2098\u2090\u2093={round(mistake(y1, y2), 2)}")
    i += 1
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.legend(framealpha=1, frameon=True)
    plt.show()
n=11
knots = np.linspace(a, b, n)
m=wyzM(knots,fun1)
derivative=[derivative_second(knots[i],fun1) for i in range(len(m))]
plt.plot(knots, derivative, color='black')
plt.scatter(knots, m, color='blue', s=40,label="m\u1D62 pochodne")
x = list(np.linspace(a, b, num=100))
i=sym.Symbol('i')
t=str(sym.diff(1/(1+i**2)))
t=str(sym.diff(t))
y = [eval(t) for i in x]
plt.plot(x,y,color='green',label="analityczna pochodna")
plt.scatter(knots, derivative, color='red',s=40,label="numeryczna pochodna")
plt.grid()
plt.title(f"Wykres 9: Drugie pochodne n={n},       h={round(knots[1]-knots[0],2)},      \u03B5\u2098\u2090\u2093={round(mistake(m,derivative),2)}")
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend(framealpha=1, frameon=True)
plt.show()