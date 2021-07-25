# -*- coding: cp1252 -*-
################################################################################
#                    M�todo da Itera��o Inversa com shift                      #
################################################################################
from numpy import *
import function_Givens

def IISK(K,mi):
    # Ordem das Matrizes
    m = len(K)
    Xkmais1 = ones([m,1],float)
    Yk = zeros([m,1],float)
    Ykmais1 = zeros([m,1],float)
    # Contador de itera��o
    n=0
    # Inicia os autovalores
    etakmais1, lbd, eta = 0, 0, 0
    #Toler�ncia de converg�ncia
    tol = 10**-6
    # Inicia a vari�vel l�gica a
    a = 0
    # C�lculo de Yk inicial
    for i in range(m):
        soma =0
        for j in range(m):
            if i==j:
                soma = soma + Xkmais1[j][0]
        Yk[i][0] = soma

    Ka = zeros([m,m],float)
    for i in range(m):
        for j in range(m):
            if i==j:
                Ka[i][j] = K[i][j] - mi
            else:
                Ka[i][j] = K[i][j]
    while a!= 1:
        if n>=1:
            eta = etakmais1
        # C�lculo de Xk+1
        Xkmais1 = function_Givens.sl(Ka,Yk)
        # C�lculo de Yk+1
        for i in range(m):
            soma =0
            for j in range(m):
                if i==j:
                    soma = soma + Xkmais1[j][0]

            Ykmais1[i][0] = soma
        # C�lculo do autovalor shifting
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

        # Contagem da itera��o
        n = n + 1
    # Atualiza eta
    eta = etakmais1

    #Retorna o autovetor e o n�mero de itera��es
    return Xkmais1,n
