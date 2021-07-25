# -*- coding: cp1252 -*-
from numpy import *
from copy import *

# Função que faz a decomposição QR.
def QRDEC(K):
    # Copia a matriz K
    R = copy(K)
    # Ordem da matriz
    n=len(R)

    # Inicia Q como sendo igual a identidade
    Q=zeros([n,n],float)
    for i in range(n):
        Q[i][i]=1

    # Fatorização
    for i in range(n):
        for h in range(i+1,n):
            if (R[i][i])**2+(R[h][i])**2==0:
                cs=1
                sn=0
            else:            
                # Determina cossenos e senos
                cs=(R[i][i])/((R[i][i])**2+(R[h][i])**2)**0.5
                sn=(R[h][i])/((R[i][i])**2+(R[h][i])**2)**0.5
            for j in range(n):
                r1,r2 = 0, 0
                r1,r2 = R[i][j], R[h][j]
                # Cálculo da matriz R
                R[i][j] = r1*cs + r2*sn
                R[h][j] = -r1*sn + r2*cs
                # Cálculo de Q
                q1,q2 = 0, 0
                q1,q2 = Q[j][i], Q[j][h]
                Q[j][i] = q1*cs + q2*sn
                Q[j][h] = -q1*sn + q2*cs
    # Retorna as matrizes Q e R
    return Q, R
