import matplotlib.pyplot as plt
import numpy as np
import random
import functools
import math
from scipy.fft import fft, ifft
def normal_distribution(x,a,c,m):
    return (a*x+c)%m
def triangular_distribution(ni,delta):
    return ni+(random.uniform(0,1)+random.uniform(0,1)-1)*delta
def distribuo(ni,delta,x):
    if x<ni:
        return -1/delta**2*(-x**2/2 +ni*x)+x/delta-(-1/delta**2*(-(ni-delta)**2/2  +ni*(ni-delta)) + (ni-delta)/delta)
    else:
        return -1/delta**2*(x**2/2 - ni*x)+x/delta-(-1/delta**2*(ni**2/2-ni**2)  +ni/delta ) + 1/2
def Chi_kwadrat(n,y,x):
    chi=0
    for i in range(len(y)):
        p=distribuo(ni,delta,y[i][0][1])-distribuo(ni,delta,y[i][0][0])
        chi+=(y[i][1]-n*p)**2/(n*p)
    return chi
x=[]
a=123
c=1
m=2**15
n=10**4
x.append(10)
for i in range(n):
    x.append(normal_distribution(x[-1],a,c,m))
x=[i/(m+1) for i in x]
plt.scatter(x[0:-1],x[1:],s=4,color="blue")
plt.xlabel(r'$X_{i}$')
plt.ylabel(r'$X_{i+1}$')
plt.title(r'U(0,1) $X_{i+1}(X_{i})$')
plt.show()
mu=sum(x)/len(x)
print(r'Średnia $ \mu $ = ',mu)
sigma1=0
for i in x:
    sigma1+=(i-mu)**2
sigma1=math.sqrt(sigma1/len(x))
print(r'Odchylenie standardowe $ \sigma $ = ',sigma1)

a = np.hstack((x))
plt.hist(a, bins=12,density=True,ec='w')
plt.xlabel(r'$X$')
plt.ylabel(r'$P(X)$')
plt.title(r'Histogram rozkładu gętości prawdopodobieństwa')
plt.show()

x=[]
a=69069
c=1
m=2**32
n=10**4
x.append(10)
for i in range(n):
    x.append(normal_distribution(x[-1],a,c,m))
x=[i/(m+1) for i in x]
plt.scatter(x[0:-1],x[1:],s=4,color="blue")
plt.xlabel(r'$X_{i}$')
plt.ylabel(r'$X_{i+1}$')
plt.title(r'U(0,1) $X_{i+1}(X_{i})$')
plt.show()

a = np.hstack((x))
plt.hist(a, bins=12,density=True,ec='w')
plt.xlabel(r'$X$')
plt.ylabel(r'$P(X)$')
plt.title(r'Histogram rozkładu gętości prawdopodobieństwa')
plt.show()

mu=sum(x)/len(x)
print(r'Średnia $ \mu $ = ',mu)
sigma1=0
for i in x:
    sigma1+=(i-mu)**2
sigma1=math.sqrt(sigma1/len(x))
print(r'Odchylenie standardowe $ \sigma $ = ',sigma1)

n=10**3
x=[]
ni=4
delta=3
K=10
for i in range(n):
    x.append(triangular_distribution(ni,delta))
p=2*delta/K
y={}
for i in range(K):
    y.setdefault(i,[(ni-delta+p*i,ni-delta+p*(i+1)),0])
for i in x:
    for j in y.values():
        if j[0][0]<=i<j[0][1]:
            j[1]+=1
print(y)
p=[]
q=[]
p2=[]
for i in range(len(y)):
    p1 = distribuo(ni, delta, y[i][0][1]) - distribuo(ni, delta, y[i][0][0])
    p2.append(y[i][1]/n)
    p.append(p1)
    q.append((y[i][0][1]+y[i][0][0])/2)
print(Chi_kwadrat(n,y,x))


u=[y[i][1] for i in range(K)]
plt.plot(range(1,K+1),u,color="black",label="f(x)")
plt.xlabel(r'$X_{i}$')
plt.ylabel(r'$X_{i+1}$')
plt.title(r'')
plt.legend(framealpha=1, frameon=True)
plt.show()
x.append(7)
x.append(1)
counts, bins = np.histogram(x)
plt.hist(bins[:-1],weights=p2,range=(1,7),ec='w')
plt.plot(q,p,'--bo',color="red")
plt.xlabel(r'$X$')
plt.ylabel(r'$n_i/N$')
plt.title(r'Histogram')
plt.show()