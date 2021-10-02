import matplotlib.pyplot as plt
import numpy as np
import copy
import random
import math
import time
from numpy.linalg import inv


def fun(x):
    return 1/(1+x**2)

def matrix_derivative(n,x):
    f=np.zeros((n+1,n+1))
    for i in range(n+1):
        f[i,0]=fun(x[i])
    for j in range(1,n+1):
        for i in range(j,n+1):
            f[i,j]=(f[i,j-1]-f[i-1,j-1])/(x[i]-x[i-j])
    return f

def polynomial_Newton(n, knots, x):
  f = matrix_derivative(n, knots)
  y = []
  for k in x:
    suma=0
    for j in range(n + 1):
      iloraz = 1
      for i in range(j):
        iloraz *= (k - knots[i])
      suma += f[j, j] * iloraz
    y.append(suma)
  return y

def Czebyszew(a,b,n):
    x=[]
    for i in range(n+1):
        x.append(1/2 * ( (a-b)*math.cos(math.pi*((2*i+1)/(2*n+2))) + (a+b)))
    return np.array(x)

n=[5,10,15,20]
colors=['y','b','c','r','g','m','tab:orange','tab:pink','tab:brown','tab:gray']
a=-5
b=5

for j in n:
    knots = np.linspace(a, b, num=j+1)
    x = list(np.linspace(a, b, num=100))
    y = polynomial_Newton(j, knots, x)
    plt.plot(x, y, color='blue',label="Interpolacja")
    plt.scatter(knots, [fun(i) for i in knots], color='red', s=40,label="węzły interpolacyjne")
    x = np.linspace(a, b, num=100)
    y = [fun(x[i]) for i in range(len(x))]
    plt.plot(x, y, color='black',label=r'f(x)=$\frac{1}{1+x^2}$')
    plt.grid()
    plt.title(f"Węzły równoodległe n={j}")
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.legend(framealpha=1, frameon=True)
    plt.show()

for j in n:
    knots = Czebyszew(a, b, j)
    x = list(np.linspace(a, b, num=100))
    y = polynomial_Newton(j, knots, x)
    plt.plot(x, y, color='blue',label="Interpolacja")
    plt.scatter(knots, [fun(i) for i in knots], color='red', s=40,label="węzły interpolacyjne")
    x = np.linspace(a, b, num=100)
    y = [fun(x[i]) for i in range(len(x))]
    plt.plot(x, y, color='black',label=r'f(x)=$\frac{1}{1+x^2}$')
    plt.grid()
    plt.title(f"Węzły Czebyszewa n={j}")
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.legend(framealpha=1, frameon=True)
    plt.show()

n=[5,6,7,8,9,10]
i=0
for j in n:
    knots = np.linspace(a, b, num=j+1)
    x = list(np.linspace(a, b, num=100))
    y = polynomial_Newton(j, knots, x)
    plt.plot(x, y, color=colors[i],label=f'n={j}')
    i+=1
plt.xlim([-5, 5])
plt.ylim([-1, 2])
x = np.linspace(a, b, num=100)
y = [fun(x[i]) for i in range(len(x))]
plt.plot(x, y, color='black',linewidth=4,label=r'f(x)=$\frac{1}{1+x^2}$')
plt.grid()
plt.title(f"Węzły równoodległe")
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend(framealpha=1, frameon=True)
plt.show()

k=0
for j in n:
    knots = Czebyszew(a, b, j)
    x = list(np.linspace(a, b, num=100))
    y = polynomial_Newton(j, knots, x)
    plt.plot(x, y, color=colors[k],label=f'n={j}')
    k+=1
plt.xlim([-5, 5])
plt.ylim([-1, 2])
x = np.linspace(a, b, num=100)
y = [fun(x[i]) for i in range(len(x))]
plt.plot(x, y, color='black',linewidth=4,label=r'f(x)=$\frac{1}{1+x^2}$')
plt.grid()
plt.title(f"Węzły Czebyszewa")
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend(framealpha=1, frameon=True)
plt.show()