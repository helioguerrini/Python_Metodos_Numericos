# -*- coding: cp1252 -*-
# Transforma uma matriz em tridiagonal
################################################################################
#                            Método de Householder                             #
################################################################################

from numpy import zeros,array,sign,identity
from copy import *
def TRDG(A):
    # Ordem da matriz
    m = len(A)
    #Inicia o vetor auxiliar
    vaux = zeros([m,1],float)
    # Inicia a matriz Ak+1
    Akmais1 = zeros([m,m],float)
    # Inicia a matriz P
    Pn = identity(m)
    for k in range(0,m-2):
        # Cálculo de alpha
        sm = 0
        for j in range(k+1,m):
            sm = sm + (A[j][k])**2
        if A[k+1][k]==0:
            alpha = -(sm**0.5)
        else:
            alpha = -sign(A[k+1][k])*(sm**0.5)
        # Cálculo de r
        r = ((0.5*alpha**2) - (0.5*alpha*A[k+1][k]))**0.5
        # Determinação do vetor w
        vetw = zeros([m-(k+1),1],float)
        vetw[0][0] = (A[k+1][k] - alpha)/(2*r)
        for j in range(k+2,m):
            vetw[j-(k+1)][0] = A[j][k]/(2*r)
        # Determina a matriz de reflexão Pn da iteração k
        for i in range(m):
            for j in range(m):
                vaux[j][0] = 0
                if j<k+1:
                    G=0
                if j==k+1:
                    G=vetw[0][0]
                if j>=k+2:
                    G=vetw[j-(k+1)][0]
                for h in range(m):
                    if h<k+1:
                        H=0
                    if h==k+1:
                        H=vetw[0][0]
                    if h>=k+2:
                        H=vetw[h-(k+1)][0]

                    if j==h:
                        vaux[j][0] = vaux[j][0] + Pn[i][h]*(1 - 2*G*H)
                    else:
                        vaux[j][0] = vaux[j][0] + Pn[i][h]*(0 - 2*G*H)
                    
            # Atualiza a linha corrente i da matriz Pn
            for l in range(m):
                Pn[i][l] = vaux[l][0]
        ### Cálculo de Ak+1 = PkAkPk ###
        # Cálculo de A'k = AkPk
        for i in range(m):
            for j in range(m):
                vaux[j][0] = 0
                if j<k+1:
                    G=0
                if j==k+1:
                    G=vetw[0][0]
                if j>=k+2:
                    G=vetw[j-(k+1)][0]
                for h in range(m):
                    if h<k+1:
                        H=0
                    if h==k+1:
                        H=vetw[0][0]
                    if h>=k+2:
                        H=vetw[h-(k+1)][0]

                    if j==h:
                        vaux[j][0] = vaux[j][0] + A[i][h]*(1 - 2*G*H)
                    else:
                        vaux[j][0] = vaux[j][0] + A[i][h]*(0 - 2*G*H)
                    
            # Atualiza a linha corrente i da matriz A
            for l in range(m):
                A[i][l] = vaux[l][0]
        # Cálculo de Ak+1 = PkA'k
        for i in range(m):
            if i<k+1:
                G=0
            if i==k+1:
                G=vetw[0][0]
            if i>=k+2:
                G=vetw[i-(k+1)][0]
                
            for j in range(m):
                Akmais1[i][j] = 0
                for h in range(m):
                    if h<k+1:
                        H=0
                    if h==k+1:
                        H=vetw[0][0]
                    if h>=k+2:
                        H=vetw[h-(k+1)][0]
                        
                    if i==h:
                        Akmais1[i][j] = Akmais1[i][j] + (1 - 2*G*H)*A[h][j]
                    else:
                        Akmais1[i][j] = Akmais1[i][j] + (0 - 2*G*H)*A[h][j]

        # Atualiza A
        A = copy(Akmais1)
    #Retorna a matriz de rigidez tridiagonal e a matriz de reflexão
    return A,Pn
