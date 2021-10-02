from math import sqrt
import numpy as np
import matplotlib.pyplot as plt

def print_matrix(a):
    for i in a:
        print (i)
    print("\n\n")

def solve_x(a,b,x):
    for i in range(5):
        w = a[i][i]
        b[i] = b[i] / w
        for j in range(5):
            a[i][j] = a[i][j] / w
        for k in range(5):
            w2 = a[k][i]
            if (k != i):
                b[k] = b[k] - b[i] * w2
                for j in range(5):
                    a[k][j] = a[k][j] - w2 * a[i][j]
    for i in range(5):
        x[i] = b[i]

def multi(a,x,c):
    for i in range(5):
        temp = 0
        for j in range(5):
            temp += a[i][j] * x[j]
        c[i] = temp;

def mean_squared_error(q):
    a = [[2 * q * pow(10, -4), 1, 6, 9, 10], [2 * pow(10, -4), 1, 6, 9, 10], [1, 6, 6, 8, 6], [5, 9, 10, 7, 10],[3, 4, 9, 7, 9]]
    a1 = [[2 * q * pow(10, -4), 1, 6, 9, 10], [2 * pow(10, -4), 1, 6, 9, 10], [1, 6, 6, 8, 6], [5, 9, 10, 7, 10],[3, 4, 9, 7, 9]]
    b = [10, 2, 9, 9, 3]
    b1 = [10, 2, 9, 9, 3]
    x = [0, 0, 0, 0, 0]
    c = [0, 0, 0, 0, 0]
    solve_x(a, b, x);
    multi(a1, x, c);
    o = 1 / 5 * sqrt(pow(c[0] - b1[0], 2) + pow(c[1] - b1[1], 2) + pow(c[2] - b1[2], 2) + pow(c[3] - b1[3], 2) + pow(c[4] - b1[4], 2))
    return o

a=[]
b=[]
for i in np.arange(0.20001, 5, 0.1):
    a.append(i)
    b.append(mean_squared_error(i))
f = open("lab1-2.txt", "w");
for i in range(len(np.arange(0.5, 5, 0.2))):
    f.write(str(a[i]))
    f.write(' ')
    f.write(str(b[i]))
    f.write('\n')
f.close()
plt.yscale('log')

plt.scatter(a,b,label="(q,o(q))",color='red')
plt.plot(a,b,label="o(q)",color="black")

plt.grid()

plt.title('Bład średniokwadratowy w zależności od q')
plt.xlabel('q')
plt.ylabel('o(q)')

plt.legend(framealpha=1, frameon=True)

plt.show()