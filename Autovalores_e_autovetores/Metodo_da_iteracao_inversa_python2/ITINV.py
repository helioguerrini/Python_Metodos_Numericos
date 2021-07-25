# -*- coding: cp1252 -*-
################################################################################
#                                                                              #
#                           Método da Iteração Inversa                         #
#                                                                              #
################################################################################
from numpy import *
import function_Givens
#### Entrada ####
# Matriz de Rigidez [K]
K = array([[11, 0, -5, -2, 0, 0],[0, 5, -4, 0, 0, 0],[-5, -4, 13, 0, -2, -2],[-2, 0, 0, 5, -3, 0],[0, 0, -2, -3, 11, -6],[0, 0, -2, 0, -6, 8]],float)
# Matriz de Massa [M]
M = array([[2, 0, 0, 0, 0, 0],[0, 5, 0, 0, 0, 0],[0, 0, 8, 0, 0, 0],[0, 0, 0, 3, 0, 0],[0, 0, 0, 0, 1, 0],[0, 0, 0, 0, 0, 2]],float)
#################
# Ordem das Matrizes
m = len(K)
Xkmais1 = ones([m,1],float)
Yk = zeros([m,1],float)
Ykmais1 = zeros([m,1],float)
# Contador de iteração
n=0
# Inicia os autovalores
lbdkmais1, lbd = 0, 0
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
while a!= 1:
    if n>=1:
        lbd = lbdkmais1
    # Cálculo de Xk+1
    Xkmais1 = function_Givens.sl(K,Yk)
    # Cálculo de Yk+1
    for i in range(m):
        soma =0
        for j in range(m):
            soma = soma + M[i][j]*Xkmais1[j][0]
        Ykmais1[i][0] = soma
    # Cálculo do autovalor
    soma1, soma2 = 0, 0
    for i in range(m):
        soma1 = soma1 + Xkmais1[i][0]*Yk[i][0]
        soma2 = soma2 + Xkmais1[i][0]*Ykmais1[i][0]
    # lbdkmais1
    lbdkmais1 = soma1/soma2
    if (((lbdkmais1 - lbd)**2)**0.5)/((lbdkmais1**2)**0.5) <= tol:
        a=1
    # Atualiza os vetores
    for i in range(m):
        Yk[i][0] = Ykmais1[i][0]
        Xkmais1[i][0] = Xkmais1[i][0]/((soma2)**0.5)
    # Contagem da iteração
    n = n + 1
# Atualiza eta
lbd = lbdkmais1
print '1° Autovalor'
print lbd
print '1° Autovetor'
print Xkmais1
print 'Número de iterações'
print n


