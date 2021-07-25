# -*- coding: cp1252 -*-
################################################################################
#          Newton’s Divided Difference Interpolating Polynomials               #
################################################################################
from numpy import array
from copy import *

def NI(listx,listFx,vpix):
    # Ordem das diferenças divididas
    m=len(list(listx))
    # Inicia a lista dos coeficientes bn para a interpolação.
    bn=copy(listFx)
    # Laço que determina a lista dos coeficientes bn para a interpolação.
    for i in range(1,m):
        for k in range(i,m): 
            bn[k]=((bn[k]-bn[i-1])/(listx[k]-listx[i-1]))
    print 'Lista dos coeficientes bn do polinômio de interpolação.'
    print bn
        
    # Determinação do valor do polinômio para o valor de x de interesse
    n=m-1
    P = bn[n]
    vpifx=copy(vpix)
    for i in range(len(vpix)):
        P = bn[n]
        x=vpix[i][0]
        for k in range(1,n+1):
            P = bn[n-k] + (x - listx[n-k])*P
        vpifx[i][0]=P

    return vpifx
