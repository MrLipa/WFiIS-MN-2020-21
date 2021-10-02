import matplotlib.pyplot as plt
import numpy as np
import copy
import random
import math
import time
from numpy.linalg import inv



def set_matrix(A):
    delta=2*L/N
    x=[0]*N
    for i in range(N):
        x[i]=-L+(i+1)*delta
        A[i,i]=pow(delta,-2)+pow(x[i],2)/2
    for i in range(1,N):
        A[i,i-1]=A[i-1,i]=-1/pow(delta,2)/2


def euclides_norm(A):
    sum=0
    for i in range(len(A)):
        sum+=A[i]*A[i]
    return math.sqrt(sum)


def QR(A):
    a=A[:,0].reshape(len(A),1)

    v=a+(euclides_norm(a)*np.eye(len(A))[:,0]).reshape(len(A),1)
    H=np.eye(len(A))-(2/(np.transpose(v)@v))*(v@np.transpose(v))
    temp=H@A
    Q=H
    R=A
    R=H@A
    for i in range(1,len(A)):
        a = temp[:, i]
        for j in range(i):
            a[j]=0
        a = a.reshape(len(A), 1)
        v=a
        if a[i]<0:
            v -= (euclides_norm(a) * np.eye(len(A))[:, i]).reshape(len(A), 1)
        else:
            v += (euclides_norm(a) * np.eye(len(A))[:, i]).reshape(len(A), 1)
        H = np.eye(len(A)) - (2 / (np.transpose(v) @ v)) * (v @ np.transpose(v))
        temp = H @ temp
        Q=Q@H
        R=H@R
    return Q,R


def eigenvalues_QR_method(A,iter=20):
    A1=copy.deepcopy(A)
    P=np.eye(len(A))
    for i in range(iter):
        Q,R=QR(A1)
        A1 = R @ Q
        P=P@Q
    val=[((np.transpose(P))@A@P)[i,i] for i in range(len(A))]
    val.sort()
    return val


def eigenvector_QR_method(A,i,iter=20):
    A1 = copy.deepcopy(A)
    P = np.eye(len(A))
    for j in range(iter):
        Q, R = QR(A1)
        A1 = R @ Q
        P = P @ Q
    H=np.transpose(P)@A@P

    x=np.full((len(A),1),0)
    i=i-1
    for j in range(len(A)-1,i,-1):
        x[j]=0
    x[i]=1
    for j in range(i-1,-1,-1):
        s=0
        for k in range(j+1,i):
            s+=H[j,k]*x[k]
        x[j]=-s/( H[j,j] - H[i,i] )
    return P@x


L=-5
N=50
iter=50

A=np.eye(N)
set_matrix(A)
colors=['k','y','b','c','r','g','m']
numbers=['\u2081','\u2082','\u2083','\u2084','\u2085']
val=eigenvalues_QR_method(A)

for j in range(5):

    x = [0] * N
    y = [0] * N
    for i in range(N):
        x[i] = i + 1
    start_time = time.time()
    y = eigenvector_QR_method(A, N-j,iter)
    print("Czas działania programu dla x=0: %s sekund" % (time.time() - start_time))
    plt.scatter(x, y, s=25, color=colors[j])
    plt.plot(x, y, label='(i,x(i))', color=colors[j])

    plt.grid()
    plt.title(u'Kolejne wartości komórek wektora własnego dla \u03BB{} = {:.2f}'.format(numbers[j], val[j]))
    # plt.xlabel('i')
    # plt.ylabel('x(i)')

    plt.legend(framealpha=1, frameon=True)

    plt.show()



start_time = time.time()
for j in range(5):
    x = [0] * N
    y = [0] * N
    for i in range(N):
        x[i] = i + 1
    y = eigenvector_QR_method(A, N-j,iter)
    plt.scatter(x, y, s=25, color=colors[j])
    plt.plot(x, y, label=u'\u03BB{} = {:.2f}'.format(numbers[j], val[j]), color=colors[j])

    plt.grid()
    plt.title("Kolejne wartości komórek 5 wektorów własnych")
    # plt.xlabel('i')
    # plt.ylabel('x(i)')

    plt.legend(framealpha=1, frameon=True)

plt.show()
print("Czas działania programu dla x=0: %s sekund" % (time.time() - start_time))
