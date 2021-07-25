# -*- coding: cp1252 -*-
# Resolução de um sistema linear ([K]{U}={R}) pelo método de fatorização de Householder
print 'Resolução de um sistema linear ([K]{U}={R}) pelo método de fatorização de Householder'
from numpy import *

# Matriz [K]
K = array([[4,-1,-1,0,-1,0],[-1,4,-1,-1,0,0],[-1,-1,4,-1,-1,-1],[0,-1,-1,4,-1,0],[-1,0,-1,-1,4,-1],[0,0,-1,0,-1,4]],float)
#K = array([[5,-4,1,0],[-4,6,-4,1],[1,-4,6,-4],[0,1,-4,5]],float)
print 'Matriz [K] = '
print K

# Vetor {b}
R=array([4,3,3,3,3,1],float)
#R=array([0,1,0,0],float)
print 'Vetor {R} = '
print R

# Ordem da matriz
n=len(K) # conta o número de colunas da matriz [K]

#inicia vetores e matrizes
P=zeros([n,n],float)
S=zeros([n,n],float)
V=zeros([n,1],float)
w=zeros([n,1],float)
ki=zeros([n,1],float)

# Fatorização de Householder
for h in range(n-1):
    norm=0
    Beta=0
    wTw=0
    # No ciclo h calcula-se o vetor Vak
    for i in range(n):
        w[i]=0
        if i>=h:
            ki[i]=K[i][h]
        else:
            ki[i]=0

        norm=norm+(ki[i])**2
    # Norma euclidiana
    norm=norm**0.5

    # Cálculo do vetor w
    for j in range(n):
        if j==h:
            w[j]=ki[j]+sign(ki[j])*norm #*e[j]=1
        else:
            w[j]=ki[j] #*e[j]=0
        wTw = wTw + w[j]**2
    #Determinação de Beta
    Beta = 2 / wTw

    # P = I - Beta x wwT do ciclo h
    for i in range(n):
        for j in range(n):
            P[i][j]=0
            P[i][j]=Beta*(P[i][j]-w[i]*w[j])
            if i==j:
                P[i][j]=P[i][j]+1

    
    # Calcula a matriz S do cilo h
    for i in range(n):
        for j in range(n):
            S[i][j]=0
            V[i]=V[i]+(P[i][j]*R[j]) # a cada ciclo h vai compondo o vetor V
            for k in range(n):
                S[i][j]=S[i][j]+(P[i][k]*K[k][j])

    # No ciclo h obtém a nova matriz [K] e o novo vetor {U}
    for i in range(n):
        R[i]=V[i]
        V[i]=0
        for j in range(n):
            K[i][j]=S[i][j]



# Resolução do sistema linear
for i in range(n):
    SS=0
    for j in range(n):
        SS=SS+(K[n-1-i][n-1-j]*V[n-1-j])
    V[n-1-i]=(R[n-1-i]-SS)/K[n-1-i][n-1-i]

print 'Solução - Vetor U = '# Solução V impresso é igual ao valor de U
print V












