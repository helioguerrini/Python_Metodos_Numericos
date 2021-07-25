# -*- coding: cp1252 -*-
################################################################################
#          Newton’s Divided Difference Interpolating Polynomials               #
################################################################################
from numpy import *
from copy import *

#################################### Entrada ###################################
# Lista de coordenadas x
listx=array([12,13,14,15,16,17,18,21,24,28,34,40,50,65,80,100])
#listx=array([-250,-200,-100,0,100,300])
# Lista de coordenadas F(x)
listFx=array([1/4.5,1/4.3,1/4.1,1/3.9,1/3.75,1/3.6,1/3.5,1/3.3,1/3.2,1/3.1,1/3.0,1/2.9,1/2.8,1/2.7,1/2.6,1/2.6])
#listFx=array([0.0163,0.318,0.699,0.87,0.941,1.04])
# Ponto de interesse
x = 101
################################################################################

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
for k in range(1,n+1):
    P = bn[n-k] + (x - listx[n-k])*P
print 'Valor do polinômio P para o valor de x ='+str(x)+' de interesse'
print P
    
