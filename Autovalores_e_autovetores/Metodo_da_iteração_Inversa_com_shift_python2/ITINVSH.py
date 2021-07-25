# -*- coding: cp1252 -*-
################################################################################
#                                                                              #
#                    Método da Iteração Inversa com shift                      #
#                                                                              #
################################################################################
from numpy import *
import function_Givens

#### Entrada ####
# Matriz de Rigidez [K]
#K = array([[11, 0, -5, -2, 0, 0],[0, 5, -4, 0, 0, 0],[-5, -4, 13, 0, -2, -2],[-2, 0, 0, 5, -3, 0],[0, 0, -2, -3, 11, -6],[0, 0, -2, 0, -6, 8]],float)
K = array([[1000,-1000,0],[-1000,2000,-1000],[0,-1000,2000]],float)

# Matriz de Massa [M]
#M = array([[2, 0, 0, 0, 0, 0],[0, 5, 0, 0, 0, 0],[0, 0, 8, 0, 0, 0],[0, 0, 0, 3, 0, 0],[0, 0, 0, 0, 1, 0],[0, 0, 0, 0, 0, 2]],float)
M = array([[0,0.6,0],[0,0.6,0],[0,0,0.6]],float)
# Valor do shifting
mi = 0.
#################
# Ordem das Matrizes
m = len(K)
Xkmais1 = ones([m,1],float)
Yk = zeros([m,1],float)
Ykmais1 = zeros([m,1],float)
# Contador de iteração
n=0
# Inicia os autovalores
etakmais1, lbd, eta = 0, 0, 0
#Tolerância de convergência
tol = 10**-6
# Inicia a variável lógica a
a = 0
# Cálculo de Yk inicial
for i in range(m):
    soma =0
    for j in range(m):
        soma = soma + M[i][j]*Xkmais1[j][0]
    Yk[i][0] = soma
Ka = zeros([m,m],float)
for i in range(m):
    for j in range(m):
        Ka[i][j] = K[i][j] - mi*M[i][j]
while a!= 1:
    if n>=1:
        eta = etakmais1
    # Cálculo de Xk+1
    Xkmais1 = function_Givens.sl(Ka,Yk)
    # Cálculo de Yk+1
    for i in range(m):
        soma =0
        for j in range(m):
            soma = soma + M[i][j]*Xkmais1[j][0]
        Ykmais1[i][0] = soma
    # Cálculo do autovalor shifting
    soma1, soma2 = 0, 0
    for i in range(m):
        soma1 = soma1 + Xkmais1[i][0]*Yk[i][0]
        soma2 = soma2 + Xkmais1[i][0]*Ykmais1[i][0]
    # etak+1
    etakmais1 = soma1/soma2
    if (((etakmais1 - eta)**2)**0.5)/((etakmais1**2)**0.5) <= tol:
        a=1
    # Atualiza os vetores
    for i in range(m):
        Yk[i][0] = Ykmais1[i][0]
        Xkmais1[i][0] = Xkmais1[i][0]/((soma2)**0.5)
    # Contagem da iteração
    n = n + 1
    print etakmais1+mi, n
# Atualiza eta
eta = etakmais1
print 'Valor de eta'
print eta
print 'Autovalor'
lbd = eta + mi
print lbd
print 'Autovetor'
print Xkmais1
print 'Número de iterações'
print n


