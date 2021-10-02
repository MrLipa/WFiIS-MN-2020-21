import matplotlib.pyplot as plt
import numpy as np
import random
import functools
import math
import scipy.special as sp
from scipy.fft import fft, ifft


def rand_normal(a,b,n):
    x=[]
    y=[]
    for i in range(n):
        u1=b-random.random()*(b-a)
        u2=b-random.random()*(b-a)
        x.append( math.sqrt(-2*math.log(u1))*math.cos(2*math.pi*u2))
        y.append( math.sqrt(-2*math.log(u1))*math.sin(2*math.pi*u2))
    return x,y


n=[10**2,10**3]
sigma=1

for h in n:
    x1,y1=rand_normal(0,1,h)
    x2,y2=rand_normal(0,1,h)
    graphx=[]
    graphy=[]
    dok=[]
    for i in np.linspace(0, 6, num=60):
        x2=[x2[j]+0.1 for j in range(h)]
        monte_carlo=0
        for k in range(h):
            monte_carlo+=1/math.sqrt(  (x1[k]-x2[k])**2  +  (y1[k]-y2[k])**2  )
        graphx.append(i)
        graphy.append(monte_carlo/h)
        dok.append(math.sqrt(math.pi)/(2*sigma)*math.exp((-i**2)/(8*sigma**2))*sp.iv(0,i ** 2 / (8 * 1 ** 2)))
    plt.plot(graphx,graphy,color="red",label=f'n={h}',linestyle='dashed')

plt.plot(graphx,dok,color="black",label=f'Wartość dokładna')
plt.title(f'Monte Carlo')
plt.legend(framealpha=1, frameon=True)
plt.show()