import matplotlib.pyplot as plt
import numpy as np
import random
import math
from scipy.fft import fft, ifft
def fun(x,m,k):
    return x**m *math.sin(k*x)
def silnia(n):
    if n>1:
        return n*silnia(n-1)
    elif n in (0,1):
        return 1;
def calka(fun,a,b,m,k,n):
    start=0
    end=0
    I=[]
    for i in range(n):
        start += (-1) ** (i) * (k * a) ** (2 * i + m + 2) / (k ** (m + 1) * silnia(2 * i + 1)*(2 * i + m + 2))
        end += (-1) ** (i) * (k * b) ** (2 * i +m + 2) / (k ** (m + 1) * silnia(2 * i + 1)*(2 * i + m + 2))
        I.append(end-start)
    return end - start,I
def calka_simpson(fun,a,b,m,k,n):
      h=(b-a)/(2*n)
      s=sum(map(lambda x: 4*fun(a+x*h,m,k) if x%2 else 2*fun(a+x*h,m,k),range(1,2*n)))+fun(a,m,k)+fun(b,m,k)
      return h/3*s
m=0
k=1
n=30
a=0
b=math.pi
C=2
c,I=calka(fun,a,b,m,k,n)
modulo=[abs(I[i]-C) for i in range(len(I))]
plt.plot(range(n),modulo,color="blue",label="|C-I|")
plt.xlabel('ilość iteracji')
plt.ylabel('|C-I|')
plt.title(f'Wykresy zależności modułu |C-I| od numeru iteracji \ndla m = 0, k = 1 (I = 2)')
plt.legend(framealpha=1, frameon=True)
plt.show()
plt.plot(range(n),I,color="black",label=r'$\int_0^{\pi} x^m sin(kx)$')
plt.xlabel('ilość iteracji')
plt.ylabel("I")
plt.title(f'Wykresy całki funkcji f(x) liczonej przy pomocy rozwinięcia w szereg \ndla m = 0, k = 1 (I = 2)')
plt.legend(framealpha=1, frameon=True)
plt.show()
m=1
k=1
C=math.pi
c,I=calka(fun,a,b,m,k,n)
modulo=[abs(I[i]-C) for i in range(len(I))]
plt.plot(range(n),modulo,color="blue",label="|C-I|")
plt.xlabel('ilość iteracji')
plt.ylabel('|C-I|')
plt.title(f'Wykresy zależności modułu |C-I| od numeru iteracji \n dla m = 1, k = 1 (I = $\pi$)')
plt.legend(framealpha=1, frameon=True)
plt.show()
plt.plot(range(n),I,color="black",label=r'$\int_0^{\pi} x^m sin(kx)$')
plt.xlabel('ilość iteracji')
plt.ylabel('I')
plt.title(f'Wykresy całki funkcji f(x) liczonej przy pomocy rozwinięcia w szereg \n dla m = 1, k = 1 (I = $\pi$)')
plt.legend(framealpha=1, frameon=True)
plt.show()
m=5
k=5
C=56.363569
c,I=calka(fun,a,b,m,k,n)
modulo=[abs(I[i]-C) for i in range(len(I))]
plt.plot(range(n),modulo,color="blue",label="|C-I|")
plt.xlabel('ilość iteracji')
plt.ylabel('|C-I|')
plt.title(f'Wykresy zależności modułu |C-I| od numeru iteracji \n dla m = 5, k = 5 (I = 56.363569)')
plt.legend(framealpha=1, frameon=True)
plt.show()
plt.plot(range(n),I,color="black",label=r'$\int_0^{\pi} x^m sin(kx)$')
plt.xlabel('ilość iteracji')
plt.ylabel('I')
plt.title(f'Wykresy całki funkcji f(x) liczonej przy pomocy rozwinięcia w szereg \n dla m = 5, k = 5 (I = 56.363569)')
plt.legend(framealpha=1, frameon=True)
plt.show()

n=[11,21,51,101,201]
m = 0
k = 1
C=2
I=[calka_simpson(fun, a, b, m, k, i) for i in n]
plt.plot(n,I,color="blue",label="|C-I|")
plt.xlabel('ilość iteracji')
plt.ylabel('|C-I|')
plt.title(f'Dokładność całki liczonej za pomoca metody Simpsona \ndla m = 0, k = 1 (I = 2)')
plt.legend(framealpha=1, frameon=True)
plt.show()
m = 1
k = 1
C=math.pi
I=[calka_simpson(fun, a, b, m, k, i) for i in n]
plt.plot(n,I,color="blue",label="|C-I|")
plt.xlabel('ilość iteracji')
plt.ylabel('|C-I|')
plt.title(f'Dokładność całki liczonej za pomoca metody Simpsona\n dla m = 1, k = 1 (I = $\pi$)')
plt.legend(framealpha=1, frameon=True)
plt.show()
m = 5
k = 5
C=56.363569
I=[calka_simpson(fun, a, b, m, k, i) for i in n]
plt.plot(n,I,color="blue",label="|C-I|")
plt.xlabel('ilość iteracji')
plt.ylabel('|C-I|')
plt.title(f'Dokładność całki liczonej za pomoca metody Simpsona\n dla m = 5, k = 5 (I = 56.363569)')
plt.legend(framealpha=1, frameon=True)
plt.show()