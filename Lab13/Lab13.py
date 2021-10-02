import matplotlib.pyplot as plt
import numpy as np
import random
import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator
def fun(x):
    return x/(4*x**2+1)
def fun1(x,k):
    return x**k
def fun2(x):
    return math.sin(x)**2
def fun3(x):
    return math.sin(x)**4
def fun4(x,y):
    np.sin(x) ** 2 * np.sin(y)**4 * math.exp(-x**2-y**2)
def fun5(x,y):
    return np.sin(x) ** 2 * np.sin(y)**4 * np.exp(-x**2-y**2)
def f(x, y):
    return np.sin(np.sqrt(x ** 2 + y ** 2))
def silnia(n):
    if n>1:
        return n*silnia(n-1)
    elif n in (0,1):
        return 1;
def calka_leandre(fun,a,b,n):
    x, w = np.polynomial.legendre.leggauss(n)
    for i in range(len(x)):
        x[i] = (b+a)/2+(b-a)/2*x[i]
    su=0
    for i in range(len(w)):
        su+=w[i]*fun(x[i])
    return su * 0.5 * (b - a), sum(w)


def calka_laguerre(fun,k,n):
    x, w = np.polynomial.laguerre.laggauss(n)
    su=0
    for i in range(len(w)):
        su+=w[i]*fun(x[i],k)
    return su,sum(w)


def calka_hermite(fun,n):
    x, w = np.polynomial.hermite.hermgauss(n)
    su=0
    for i in range(len(w)):
        su+=w[i]*fun(x[i])
    return su

a=2
c=1
f=[]
w=[]
c1=1/(2*2**2)*math.log(2**2*2**2+c**2)-1/(2*2**2)*math.log(2**2*0**2+c**2)
for i in range(2,20):
    f.append(abs(calka_leandre(fun,0,2,i)[0]-c1))
    w.append(calka_leandre(fun,0,2,i)[1])
plt.plot(range(2,20),f,color="black",label=r'$|c_1-c_{1,a}|$')
plt.xlabel('n')
plt.ylabel(r'$|c_1-c_{1,a}|$')
plt.title(r'Wykresy zależności n od $|c_1-c_{1,a}|$')
plt.legend(framealpha=1, frameon=True)
plt.show()

x=np.linspace(0, 2, num=100)
y=[fun(i) for i in x]
plt.plot(x,y,color="black",label="f(x)")
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title(r'Wykres funkcji $f(x)=\frac{x}{4x^2+1}$')
plt.legend(framealpha=1, frameon=True)
plt.show()

f=[]
w=[]
k=5
c2=silnia(k)
for i in range(2,20):
    f.append(abs(calka_laguerre(fun1,k,i)[0]-c2))
    w.append(calka_laguerre(fun1,k,i)[1])
plt.plot(range(2,20),f,color="black",label=r'$|c_2-c_{2,a}|$')
plt.xlabel('n')
plt.ylabel(r'$|c_2-c_{2,a}|$')
plt.title(r'Wykresy zależności n od $|c_2-c_{2,a}|$ dla k=5')
plt.legend(framealpha=1, frameon=True)
plt.show()

f=[]
w=[]
k=10
c2=silnia(k)
for i in range(2,20):
    f.append(abs(calka_laguerre(fun1,k,i)[0]-c2))
    w.append(calka_laguerre(fun1,k,i)[1])
plt.plot(range(2,20),f,color="black",label=r'$|c_2-c_{2,a}|$')
plt.xlabel('n')
plt.ylabel(r'$|c_2-c_{2,a}|$')
plt.title(r'Wykresy zależności n od $|c_2-c_{2,a}|$ dla k=10')
plt.legend(framealpha=1, frameon=True)
plt.show()

x=np.linspace(0, 2, num=100)
y=[fun1(i,k) for i in x]
plt.plot(x,y,color="black",label="f(x)")
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title(r'Wykres funkcji $f(x)=x^kexp(-x)$')
plt.legend(framealpha=1, frameon=True)
plt.show()

f=[]
c3=0.1919832644
for i in range(2,15):
    f.append(abs((calka_hermite(fun2,i)*calka_hermite(fun3,i))-c3))
plt.plot(range(2,15),f,color="black",label=r'$|c_3-c_{dok}|$')
plt.xlabel('n')
plt.ylabel(r'$|c_3-c_{dok}|$')
plt.title(r'Wykresy zależności n od $|c_3-c_{dok}|$')
plt.legend(framealpha=1, frameon=True)
plt.show()
x=np.linspace(-5, 5, num=100)
y=np.linspace(-5, 5, num=100)
z=[fun4(i[0],i[1]) for i in zip(x,y)]

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

x = np.arange(-5, 5, 0.25)
y = np.arange(-5, 5, 0.25)
x, y = np.meshgrid(x, y)
z = fun5(x, y)

surf = ax.plot_surface(x, y, z, cmap=cm.coolwarm,linewidth=0, antialiased=False)

ax.set_zlim(0.0, 0.08)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter('{x:.02f}')

fig.colorbar(surf, shrink=0.5, aspect=5)
plt.title(r'Wykres funkcji $f(x,y)=sin^2(x)cos^4(x)exp(-x^2-y^2)$')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("f(x,y)")

plt.show()