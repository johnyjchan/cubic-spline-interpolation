# -*- coding:utf-8 -*-
'''
Created on 2017/04/28

@author: John Chan
@version: 1.0
@contact: chenyj@bcc.ac.cn

Cubic Spline Interpolation provides numeric computing formula to interpolate curve.
This source code was designed to draw a 3D curve.
If you have any question or optimized idea, welcome to contact me.

'''
import math
import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from numpy import matrix, average
import scipy.linalg 

# Parameters
pointsInterpolation=False
curveInterpolation=True
'''
    numberOfInterpolation determines the precision of interpolation.
    bigger numberOfInterpolation, more smooth curve
'''
numberOfInterpolation = 100


j=0
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
       
def cubicSplineInterpolate(x_axis,y_axis,z_axis):
    dx=[]
    dy=[]
    dz=[]
    matrix=[]
    n=2
    while n<len(x_axis):
        dx.append(3*(x_axis[n]-2*x_axis[n-1]+x_axis[n-2]))
        dy.append(3*(y_axis[n]-2*y_axis[n-1]+y_axis[n-2]))
        dz.append(3*(z_axis[n]-2*z_axis[n-1]+z_axis[n-2]))
        n=n+1
        
    '''
        produce square matrix looks like :
        [[2.0, 0.5, 0.0, 0.0], [0.5, 2.0, 0.5, 0.0], [0.0, 0.5, 2.0, 0.5], [0.0, 0.0, 2.0, 0.5]]
        the classes of the matrix depends on the length of x_axis(number of nodes)
    '''
    matrix.append([float(2), float(0.5)])
    for m in range(len(x_axis)-4):
        matrix[0].append(float(0))                
    n=2
    while n<len(x_axis)-2:
        matrix.append([])
        for m in range(n-2):
            matrix[n-1].append(float(0)) 
              
        matrix[n-1].append(float(0.5))
        matrix[n-1].append(float(2))
        matrix[n-1].append(float(0.5))
        
        for m in range(len(x_axis)-n-3):
            matrix[n-1].append(float(0)) 
        n=n+1
        
    matrix.append([])
    for m in range(n-2):
        matrix[n-1].append(float(0))    
    matrix[n-1].append(float(0.5))    
    matrix[n-1].append(float(2))
    

    '''
        LU Factorization may not be optimal method to solve this regular matrix.
        If you guys have better idea to solve the Equation, please contact me.
    '''
    P, L, U = doLUFactorization(matrix)
    u=solveEquations(P,L,U,dx)
    v=solveEquations(P,L,U,dy)
    w=solveEquations(P,L,U,dz)
    
    '''
        prepare right-side vector
    '''
    m=0
    U=[0]
    V=[0]
    W=[0]
    while m<len(u):
        U.append(u[m])
        V.append(v[m])
        W.append(w[m])
        m=m+1
    U.append(0)
    V.append(0)
    W.append(0)
   
    plotCubicSpline(U,V,W,x_axis,y_axis,z_axis)

'''
    calculate each parameters of location.
'''
def func(x1,x2,t,v1,v2,t1,t2):
    ft=((t2-t)**3*v1+(t-t1)**3*v2)/6+(t-t1)*(x2-v2/6)+(t2-t)*(x1-v1/6)
    return ft

'''
    note: 
    too many points interpolate points make your computer slack.
    To interpolate large number of input parameters,
    please switch to ax.plot().
'''
def plotCubicSpline(U,V,W,x_axis,y_axis,z_axis):
    m=1
    xLinespace=[]
    yLinespace=[]
    zLinespace=[]
    while m<len(x_axis):
        for t in np.arange(m-1,m,1/float(numberOfInterpolation)):
            xLinespace.append(func(x_axis[m-1],x_axis[m],t,U[m-1],U[m],m-1,m))
            yLinespace.append(func(y_axis[m-1],y_axis[m],t,V[m-1],V[m],m-1,m))
            zLinespace.append(func(z_axis[m-1],z_axis[m],t,W[m-1],W[m],m-1,m))
        m=m+1
    if pointsInterpolation:
        ax.scatter(xLinespace, yLinespace,zLinespace,color="red",s=0.01)
    if curveInterpolation:
        ax.plot(xLinespace, yLinespace,zLinespace,color="red")
    '''
    matched group, annotate it if unnecessary
    '''
    ax.plot(x_axis,y_axis,z_axis,color="blue")

'''
    matrix·x = y
    P·matrix = L·U
    P·matrix·x = L·U·x = P·y
    L·U·x = y1
    U·x = y2
    x = y3
'''  
def solveEquations(P,L,U,y):
    y1=np.dot(P,y)
    y2=y1
    m=0
    while m<len(y):
        n=0
        while n<m:
            y2[m]=y2[m]-y2[n]*L[m][n]
            n=n+1
        y2[m]=y2[m]/L[m][m]
        m=m+1
    y3=y2
    m=len(y)-1
    while m>=0:
        n=len(y)-1
        while n>m:
            y3[m]=y3[m]-y3[n]*U[m][n]
            n=n-1
        y3[m]=y3[m]/U[m][m]
        m=m-1
    return y3

'''
    this is the Scipy tool with high complexity.
'''    
def doLUFactorization(matrix):    
    P, L, U=scipy.linalg.lu(matrix)
    return P, L, U   

'''
    input parameters
    each linespace contians at least 3 elements
'''
x_axis = [1, 2, 3, 4]
y_axis = [2, 3, 4, 5]
z_axis = [3, 4, 7, 5]

cubicSplineInterpolate(x_axis,y_axis,z_axis)

plt.show()
