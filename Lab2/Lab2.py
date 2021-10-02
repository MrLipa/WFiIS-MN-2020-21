import matplotlib.pyplot as plt
import numpy as np
import copy
import random

def print_matrix(A,nazwa,rount=1):
    if type(A[0])==int or type(A[0])==float:
        print(f"__________Wektor {nazwa}[]__________")
        for i in range(N):
            print(round(A[i],rount))
        print()
        print()
    if type(A[0])==list:
        print(f"__________Macierz {nazwa}[]__________")
        for i in range(N):
            round_to_whole = [round(num, rount) for num in A[i]]
            print(f'[{round_to_whole}]')
        print()
        print()

def print_equation_matrix(A,x,b,rount=1):
    print(f"__________A*x=b__________")
    print()
    temp=int(N/2)
    for i in range(N):
        round_to_whole = [round(num, rount) for num in A[i]]
        if i==temp:
            print(round_to_whole, '* [', round(x[i],rount), ']', " = ", '[', round(b[i],rount), ']')
        else:
            print(round_to_whole, '  [', round(x[i],rount), ']', "   ", '[', round(b[i],rount), ']')
    print()
    print()



def set_matrix(A,c,key=0):
    if key==0:
        for i in range(N):
            for j in range(N):
                if j==0:
                    A[i][j]=1
                elif j==1:
                    A[i][j]=random.uniform(-1.,1.)
                else:
                    A[i][j] = A[i][1] ** j
        for i in range(N):
            c[i]=random.randint(-10, 10)
        if np.linalg.det(A) == 0:
            set_matrix(A, c)
        else:
            return
    elif key==1:
        N1 = 300
        y1 = [0 for i in range(N1)]
        x1 = [0 for i in range(N1)]
        for i in range(0, N1):
            x1[i] =random.uniform(-1.,1.)
            for j in range(N):
                y1[i] += ((x1[i]) ** j) * c1[j]
        return x1,y1

def LU(A):
    L = [[0 for i in range(N)] for i in range(N)]
    for i in range(0, N):
        L[i][i] = 1
    U = [[0 for i in range(N)] for i in range(N)]

    for j in range(N):
        for i in range(j + 1):
            s = 0
            for k in range(i):
                s += A[i][k] * A[k][j]
            A[i][j] -= s
            U[i][j]=A[i][j]
        for i in range(j + 1, N):
            s = 0
            for k in range(j):
                s += A[i][k] * A[k][j]
            A[i][j] = (A[i][j] - s) / U[j][j]
            L[i][j] = A[i][j]
    return (L,U)


def calculate_LU(A,x,b):
    A1=copy.deepcopy(A)
    L, U = LU(A1)
    z=[]
    for i in range(N):
        z.append(0)
    for i in range(N):
        temp = 0
        for j in range(N):
            temp+=z[j]*A1[i][j]
        z[i]=b[i]-temp
    for i in range(N-1,-1,-1):
        temp=0
        for j in range(N-1,-1,-1):
            temp+=x[j]*A1[i][j]

        x[i]=(z[i]-temp)/U[i][i]


def multiply_matrix(A,x,b):
    for i in range(N):
        for j in range(N):
            b[i] += A[i][j] * x[j]


def horner(x,wsp):
    wynik = wsp[2]
    for i in range(N-1):
        wynik = wynik * x + wsp[N-i-2]

    return wynik

def horner_matrix(x,c,b):
    for i in range(N):
        b[i]=horner(x[i][1],c)


def determinant_matrix_LU(A):
    A1 = copy.deepcopy(A)
    L, U = LU(A1)
    temp=1
    for i in range(N):
        temp*=U[i][i]
    return temp

def transposition_matrix(A):
    B=[]
    for i in range(N):
        B.append([0])
        for j in range(N - 1):
            temp = B[i]
            temp.append(0)
    for iw in range(N):
        for ik in range(N):
            B[ik][iw] = A[iw][ik]
    return B


def inverse_matrix(A):
    A1 = copy.deepcopy(A)
    L, U = LU(A1)
    A2=[]
    for i in range(N):
        A2.append([0])
        for j in range(N - 1):
            temp = A2[i]
            temp.append(0)
    for i in range(N):
        B=[0]*N
        I=[0]*N
        I[i]=1
        calculate_LU(A,B,I)
        A2[i]=B
    return transposition_matrix(A2)


def index_matrix_conditioning(A):
    B=[]
    for i in range(N):
        temp=0
        for j in range(N):
            temp+=abs(A[j][i])
        B.append(temp)
    A1=inverse_matrix(A)
    B1 = []
    for i in range(N):
        temp = 0
        for j in range(N):
            temp += abs(A1[j][i])
        B1.append(temp)
    return max(B)*max(B1)

N = 20
A=[0]*N
for i in range(N):
    A[i]=[0]*N
c=[0]*N
y=[0]*N

set_matrix(A,c)
print_equation_matrix(A,c,y)
horner_matrix(A,c,y)
print_equation_matrix(A,c,y)

A1=copy.deepcopy(A)
c1=[0]*N
print_equation_matrix(A,c1,y)
calculate_LU(A1,c1,y)
print_equation_matrix(A,c1,y)

x1,y1=set_matrix(A,c1,1)
x=[A[i][1] for i in range(N)]

plt.scatter(x1,y1,label="(x1,y1)",s=4,color='blue')
plt.scatter(x,y,label="(x,y)",s=25,color='red')

plt.grid()

plt.xlabel('x')
plt.ylabel('y(x)')

plt.legend(framealpha=1, frameon=True);

plt.show()


print_matrix(inverse_matrix(A),'A\'')
print("Wyznacznik macierzy A: ",determinant_matrix_LU(A))
print("Współczynnik uwarunkowania macierzy A: ",round(index_matrix_conditioning(A),2))
