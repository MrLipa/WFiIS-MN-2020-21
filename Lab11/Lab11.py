import matplotlib.pyplot as plt
import numpy as np
import random
import math
from scipy.fft import fft, ifft

def fun1(t,omega):
    return math.sin(1*omega*t) + math.sin(2*omega*t) + math.sin(3*omega*t)
def fun2(t,sigma):
    return 1/(sigma*math.sqrt(2*math.pi))*math.exp(-t**2/(2*sigma**2))
def Fourier(k,T,a,b,n):
    N=2**k
    Tmax=3*T
    sigma=T/n
    omega = 2 * math.pi / T
    dt = Tmax / N
    t = np.arange(0, Tmax, dt)
    f0=[fun1(i,omega) for i in t]
    f = [fun1(i, omega)+random.uniform(-0.5,0.5) for i in t]
    g1 = [fun2(i, sigma) for i in t]
    g2 = [fun2(i, sigma) for i in t]
    g=[g1[i]+g2[i] for i in range(len(g1))]
    FFT_F = fft(f)
    FFT_G = fft(g)
    f1=ifft([FFT_F[i]*FFT_G[i] for i in range(len(FFT_F))])
    Fmax = max([abs(i) for i in f1])

    plt.plot(t,f0,color="black",label="f\u2080(t\u1D62)")
    plt.plot(t, f, color="grey", label="f(t\u1D62)=f\u2080(t\u1D62)+\u0394",linewidth=0.3)
    plt.plot(t, f1 * 2.5 / Fmax, color='red', label="f*g")
    plt.xlabel('t')
    plt.ylabel('f(t)')
    plt.title(f'Wykresy funkcji sygnału niezaburzonego,\nzaburzonego i znormalizowanego splotu dla k={k}')
    plt.legend(framealpha=1, frameon=True)
    plt.show()
k=[8,10,12]
T=1
a=-1.5
b=1.5
n=20
for i in k:
    Fourier(i,T,a,b,n)
# n=10
# for i in k:
#     Fourier(i,T,a,b,n)
# n=100
# for i in k:
#     Fourier(i,T,a,b,n)