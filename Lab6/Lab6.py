import matplotlib.pyplot as plt
import numpy as np
import copy
import random
import math
import time
from numpy.linalg import inv

def delta_r(x,y):
    r=np.zeros((2,2))
    temp=np.zeros((2,1))
    temp[0,0]=2*x*y**2-3*x**2*y-2
    temp[1,0]=x**2*y**3+2*x*y-12
    r[0,0]=2*y**2 -6*x*y
    r[1,0]=2*x*y**3+2*y
    r[0,1]=4*x*y-3*x**2
    r[1,1]=3*x**2*y**2+2*x
    r=-np.linalg.inv(r)@temp
    return r

# def derivative(x,delta=0.001):
#     return (fun(x+delta)-fun(x - delta))/(2 * delta)

def euclides_norm(A):
    sum=0
    for i in range(len(A)):
        sum+=A[i,0]*A[i,0]
    return math.sqrt(sum)

def Newton_metod(r):
    k=[]
    delta_x=[]
    i=0
    x=[]
    y=[]
    x.append(r[0,0])
    y.append(r[1,0])
    while True :
        i+=1
        k.append(i)
        temp=r
        delta_x.append(euclides_norm(delta_r(r[0,0],r[1,0])))
        r=r+delta_r(r[0,0],r[1,0])
        if euclides_norm(r-temp)<10e-6:
            break
        x.append(r[0,0])
        y.append(r[1,0])
    return r,k,delta_x,x,y

r=np.zeros((2,1))
r[0,0]=10
r[1,0]=10
r,k,delta_x,x,y=Newton_metod(r)

plt.yscale('log')
plt.plot(k,delta_x,color='blue')
plt.scatter(k,delta_x,s=15,color='black')
plt.grid()
plt.xlabel('k')
plt.ylabel('\u0394 x')
plt.show()
print("Ilość iteracji",len(k))
print("Rozwiązanie",r[0,0],r[1,0])
for i in range(len(x)-1):
    plt.arrow(x[i],y[i],-math.sqrt((x[i]-x[i+1])**2),-math.sqrt((y[i]-y[i+1])**2),length_includes_head=True,head_width=0.3, head_length=0.4)
plt.annotate(f"[{round(r[0,0])},{round(r[1,0])}]", (1, 2),xytext =(3, 2), arrowprops = dict(facecolor ='green',shrink = 0.02))
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.show()

r=np.zeros((2,1))
r[0,0]=10
r[1,0]=-4
r,k,delta_x,x,y=Newton_metod(r)

plt.yscale('log')
plt.plot(k,delta_x,color='blue')
plt.scatter(k,delta_x,s=15,color='black')
plt.grid()
plt.xlabel('k')
plt.ylabel('\u0394 x')
plt.show()
plt.scatter(x,y,s=15,color='blue')
print("Ilość iteracji",len(k))
print("Rozwiązanie",r[0,0],r[1,0])
plt.annotate(f"[{round(r[0,0])},{round(r[1,0])}]", (1, 2),xytext =(4, 1), arrowprops = dict(facecolor ='green',shrink = 0.02))
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
for i in range(len(x)-1):
    if x[i] - x[i + 1] > 0 and y[i] - y[i + 1] > 0:
        plt.arrow(x[i],y[i],-math.sqrt((x[i]-x[i+1])**2),-math.sqrt((y[i]-y[i+1])**2),length_includes_head=True,head_width=0.3, head_length=0.4)
    if x[i] - x[i + 1] < 0 and y[i] - y[i + 1] < 0:
        plt.arrow(x[i],y[i],math.sqrt((x[i]-x[i+1])**2),math.sqrt((y[i]-y[i+1])**2),length_includes_head=True,head_width=0.3, head_length=0.4)
    if x[i] - x[i + 1] > 0 and y[i] - y[i + 1] < 0:
        plt.arrow(x[i],y[i],-math.sqrt((x[i]-x[i+1])**2),math.sqrt((y[i]-y[i+1])**2),length_includes_head=True,head_width=0.3, head_length=0.4)
    if x[i] - x[i + 1] < 0 and y[i] - y[i + 1] > 0:
        plt.arrow(x[i],y[i],math.sqrt((x[i]-x[i+1])**2),-math.sqrt((y[i]-y[i+1])**2),length_includes_head=True,head_width=0.3, head_length=0.4)
plt.show()

r=np.zeros((2,1))
r[0,0]=50
r[1,0]=50
r,k,delta_x,x,y=Newton_metod(r)

plt.scatter(x,y,s=15,color='blue')
print("Ilość iteracji",len(k))
print("Rozwiązanie",r[0,0],r[1,0])
plt.annotate(f"[{round(r[0,0])},{round(r[1,0])}]", (1, 2),xytext =(10, 5), arrowprops = dict(facecolor ='green',shrink = 0.02))
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
for i in range(len(x)-1):
    if x[i] - x[i + 1] > 0 and y[i] - y[i + 1] > 0:
        plt.arrow(x[i],y[i],-math.sqrt((x[i]-x[i+1])**2),-math.sqrt((y[i]-y[i+1])**2),length_includes_head=True,head_width=1.2, head_length=1)
    if x[i] - x[i + 1] < 0 and y[i] - y[i + 1] < 0:
        plt.arrow(x[i],y[i],math.sqrt((x[i]-x[i+1])**2),math.sqrt((y[i]-y[i+1])**2),length_includes_head=True,head_width=1.2, head_length=1)
    if x[i] - x[i + 1] > 0 and y[i] - y[i + 1] < 0:
        plt.arrow(x[i],y[i],-math.sqrt((x[i]-x[i+1])**2),math.sqrt((y[i]-y[i+1])**2),length_includes_head=True,head_width=1.2, head_length=1)
    if x[i] - x[i + 1] < 0 and y[i] - y[i + 1] > 0:
        plt.arrow(x[i],y[i],math.sqrt((x[i]-x[i+1])**2),-math.sqrt((y[i]-y[i+1])**2),length_includes_head=True,head_width=1.2, head_length=1)
plt.show()

r=np.zeros((2,1))
r[0,0]=-30
r[1,0]=-30
r,k,delta_x,x,y=Newton_metod(r)

plt.scatter(x,y,s=15,color='blue')
print("Ilość iteracji",len(k))
print("Rozwiązanie",r[0,0],r[1,0])
plt.annotate(f"[{round(r[0,0])},{round(r[1,0])}]", (1, 2),xytext =(-45, -260), arrowprops = dict(facecolor ='green',shrink = 0.02))
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
for i in range(len(x)-1):
    if x[i] - x[i + 1] > 0 and y[i] - y[i + 1] > 0:
        plt.arrow(x[i],y[i],-math.sqrt((x[i]-x[i+1])**2),-math.sqrt((y[i]-y[i+1])**2),length_includes_head=True,head_width=10, head_length=12)
    if x[i] - x[i + 1] < 0 and y[i] - y[i + 1] < 0:
        plt.arrow(x[i],y[i],math.sqrt((x[i]-x[i+1])**2),math.sqrt((y[i]-y[i+1])**2),length_includes_head=True,head_width=10, head_length=12)
    if x[i] - x[i + 1] > 0 and y[i] - y[i + 1] < 0:
        plt.arrow(x[i],y[i],-math.sqrt((x[i]-x[i+1])**2),math.sqrt((y[i]-y[i+1])**2),length_includes_head=True,head_width=10, head_length=12)
    if x[i] - x[i + 1] < 0 and y[i] - y[i + 1] > 0:
        plt.arrow(x[i],y[i],math.sqrt((x[i]-x[i+1])**2),-math.sqrt((y[i]-y[i+1])**2),length_includes_head=True,head_width=10, head_length=12)
plt.show()