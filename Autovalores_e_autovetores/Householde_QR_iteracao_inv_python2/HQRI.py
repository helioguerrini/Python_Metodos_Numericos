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
#    Transforma a matriz K numa tridiagonal pelo m�todo de Householder      #
#############################################################################
KPn = function_HSEHLD2.TRDG(K)
T, Pn = KPn[0], KPn[1]
print 'Matriz K tridiagonal'
print T
print 'Matriz Pn de reflex�o'
print Pn
# Armazena a matriz tridiagonal
K = copy(T)
KPn=[]
#############################################################################
#                           M�todo da itera��o QR                           #
#############################################################################
# Ordem da matriz
m = len(T)
# Contador de itera��o e vari�vel l�gica
I, a = 0, 0
# Toler�ncia a converg�ncia
tol = 10**-9
# Inicia as normas de autovalores para o controle de converg�ncia
lbd, lbdImais1 = 0, 0
########################### La�o da itera��o QR ############################
while a!= 1:
    if I>=1:
        lbd = lbdImais1
    # Obt�m as matrizes Q e R
    MQR = function_QR.QRDEC(T)
    Q, R = MQR[0], MQR[1]
    # C�lculo de Tmais1 
    for i in range(m):
        for j in range(m):
            T[i][j]=0
            for h in range(m):
                # A matriz T � atualizada para Tmais1 
                T[i][j] = T[i][j] + R[i][h]*Q[h][j]
    # Para o controle da converg�ncia toma-se o 1� autovalor
    lbdImais1 = T[m-1][m-1]
    # Verifica-se a converg�ncia
    if (((lbdImais1 - lbd)**2)**0.5)/((lbdImais1**2)**0.5) <= tol or I>500:
        a=1
    I=I+1

# Imprime resultados
print 'N�mero de itera��es do m�todo QR =', I
print 'Matriz T que cont�m os autovalores na diagonal principal ='
print T

# Armazena os autovalores na lista de autovalores
ListEgvl=[]
for i in range(m):
    ListEgvl.append(T[(m-1)-i][(m-1)-i])
#############################################################################
#                   M�todo da Itera��o Inversa com shift                    #
#############################################################################
# Com o m�todo da itera��o iversa acha-se os autovetores da matriz transformada
#em tridiagonal. Multiplicando cada autovetor encontrado pela matriz de reflex�o
#Pn encontra-se o respectivo autovetor da matriz K original.
# Cada autovalor encontrado anteriormente � utilizado como shift.

#Inicia o autovetor
eigvtr = zeros([m,1],float)
for i in range(m):
    mi = ListEgvl[i]
    KPn = function_ITIS.IISK(K,mi)
    # Produto entre a Matriz de reflex�o Pn e o autovetor de T
    # Tal produto serve para se encontrar o autovetor de K
    for j in range(m):
        eigvtr[j][0] = 0
        for h in range(m):
            eigvtr[j][0] = eigvtr[j][0] + Pn[j][h]*KPn[0][h][0]
    print 'N�mero de itera��es inversas =',KPn[1]
    print 'Autovalor =', ListEgvl[i]
    print 'Autovetor ='
    print eigvtr

