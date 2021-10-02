import matplotlib.pyplot as plt
import numpy as np

N=100

def wypisz(a):
    for i in range(N):
        print(a[i])
    print()
    print()


def ustaw(a, omega, h):
    a[0][0] = 1
    a[1][0] = -1
    for i in range(N):
        a[i][i] = 1
        if i >= 2:
            a[i][i - 2] = 1
            a[i][i - 1] = omega * omega * h * h - 2


def policz(x, b, a):
    for i in range(N):
        w = a[i][i]
        b[i] = b[i] / w
        for j in range(N):
            a[i][j] = a[i][j] / w
        for k in range(N):
            w2 = a[k][i]
            if k != i:
                b[k] = b[k] - b[i] * w2
                for j in range(N):
                    a[k][j] = a[k][j] - w2 * a[i][j]
    for i in range(N):
        x[i] = b[i]

A = 1
v_0 = 0
omega = 1
h = 0.1

a = []
b = []
x = []
y = []

for i in range(N):
    b.append(0)
for i in range(N):
    x.append(0)
for i in range(N):
    y.append(0)
for i in range(N):
    a.append([0])
    for j in range(N-1):
        temp=a[i]
        temp.append(0)

b[0]=A
b[1]=v_0 * h

for i in range(N):
    y[i]=h*i

x1=list(np.arange(0,N*0.1,0.1))
y1 =list(np.cos(x1))


ustaw(a,omega,h)
policz(x,b,a)

f = open("lab1-1.txt", "w");
for i in range(N):
    f.write(str(y[i]))
    f.write(' ')
    f.write(str(x[i]))
    f.write('\n')
f.close()

plt.plot(y, x,label="x(t)",color="blue")
plt.plot(x1,y1,label="cos(t)",color="black")
plt.scatter(y,x,label="dt = 0.1, x(t)",color='red')

plt.title('Wychylenie x(t)')
plt.xlabel('t')
plt.ylabel('x(t)')

plt.grid()

plt.legend(framealpha=1, frameon=True);

plt.show()


