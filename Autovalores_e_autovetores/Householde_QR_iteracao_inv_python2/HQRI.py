# -*- coding: cp1252 -*-
#############################################################################
#               Householder-QR-Inverse Iteration (HQRI) solution            #
#############################################################################
from numpy import zeros,array
from copy import *
import function_QR
import function_HSEHLD2
import function_ITIS
################################ ENTRADA ####################################
# Matriz [K]
#K = array([[5.5,0,-1.25,-0.816496581,0,0],[0,1,-0.632455532,0,0,0],[-1.25,-0.632455532,1.625,0,-0.707106781,-0.5],[-0.816496581,0,0,1.666666667,-1.732050808,0],[0,0,-0.707106781,-1.732050808,11,-4.242640687],[0,0,-0.5,0,-4.242640687,4]],float)
K = array([[1666.6666666,-1666.6666666,0],[-1666.6666666,3333.3333333,-1666.6666666],[0,-1666.6666666,3333.333333333]],float)
#############################################################################
#    Transforma a matriz K numa tridiagonal pelo método de Householder      #
#############################################################################
KPn = function_HSEHLD2.TRDG(K)
T, Pn = KPn[0], KPn[1]
print 'Matriz K tridiagonal'
print T
print 'Matriz Pn de reflexão'
print Pn
# Armazena a matriz tridiagonal
K = copy(T)
KPn=[]
#############################################################################
#                           Método da iteração QR                           #
#############################################################################
# Ordem da matriz
m = len(T)
# Contador de iteração e variável lógica
I, a = 0, 0
# Tolerância a convergência
tol = 10**-9
# Inicia as normas de autovalores para o controle de convergência
lbd, lbdImais1 = 0, 0
########################### Laço da iteração QR ############################
while a!= 1:
    if I>=1:
        lbd = lbdImais1
    # Obtém as matrizes Q e R
    MQR = function_QR.QRDEC(T)
    Q, R = MQR[0], MQR[1]
    # Cálculo de Tmais1 
    for i in range(m):
        for j in range(m):
            T[i][j]=0
            for h in range(m):
                # A matriz T é atualizada para Tmais1 
                T[i][j] = T[i][j] + R[i][h]*Q[h][j]
    # Para o controle da convergência toma-se o 1° autovalor
    lbdImais1 = T[m-1][m-1]
    # Verifica-se a convergência
    if (((lbdImais1 - lbd)**2)**0.5)/((lbdImais1**2)**0.5) <= tol or I>500:
        a=1
    I=I+1

# Imprime resultados
print 'Número de iterações do método QR =', I
print 'Matriz T que contém os autovalores na diagonal principal ='
print T

# Armazena os autovalores na lista de autovalores
ListEgvl=[]
for i in range(m):
    ListEgvl.append(T[(m-1)-i][(m-1)-i])
#############################################################################
#                   Método da Iteração Inversa com shift                    #
#############################################################################
# Com o método da iteração iversa acha-se os autovetores da matriz transformada
#em tridiagonal. Multiplicando cada autovetor encontrado pela matriz de reflexão
#Pn encontra-se o respectivo autovetor da matriz K original.
# Cada autovalor encontrado anteriormente é utilizado como shift.

#Inicia o autovetor
eigvtr = zeros([m,1],float)
for i in range(m):
    mi = ListEgvl[i]
    KPn = function_ITIS.IISK(K,mi)
    # Produto entre a Matriz de reflexão Pn e o autovetor de T
    # Tal produto serve para se encontrar o autovetor de K
    for j in range(m):
        eigvtr[j][0] = 0
        for h in range(m):
            eigvtr[j][0] = eigvtr[j][0] + Pn[j][h]*KPn[0][h][0]
    print 'Número de iterações inversas =',KPn[1]
    print 'Autovalor =', ListEgvl[i]
    print 'Autovetor ='
    print eigvtr

