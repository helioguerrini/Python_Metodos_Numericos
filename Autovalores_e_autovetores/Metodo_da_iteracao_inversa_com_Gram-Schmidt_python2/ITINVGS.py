# -*- coding: cp1252 -*-
################################################################################
#                                                                              #
#       M�todo da Itera��o Inversa com Ortonormaliza��o de Gram-Schmidt        #
#                                                                              #
################################################################################
from numpy import *
import function_Givens
import function_ITINV

#### Entrada ####

# Matriz de Rigidez [K]
#K = array([[79.863, 2.653, 18.568],[2.653, 0, 0],[18.568, 0, 106.67]],float)
K = array([[1000, -1000, 0],[-1000, 2000, -1000],[0, -1000, 2000]],float)
# Matriz de Massa [M]
#M = array([[1, 0, 0],[0, 1, 0],[0, 0, 1]],float)
M = array([[0.6, 0, 0],[0,0.6, 0],[0, 0, 0.6]],float)
#Toler�ncia de converg�ncia
tol = 10**-6
# N�mero de autopares requisitados
ap = 3
#################
# Ordem das Matrizes
m = len(K)
# Inicia vetores
Xkmais1 = ones([m,1],float)
Xotkmais1 = ones([m,1],float)
Yk = zeros([m,1],float)
# Lista de autovetores
listegvectors = []
# Resultados da itera��o inversa
RESINTINV = function_ITINV.eigenvv(K,M,Xkmais1,m,tol)
# Adiciona o primeiro autovetor na lista de autovetores
listegvectors.append(RESINTINV[1])
# Imprime resultados da itera��o inversa
print '1� Autovalor =',RESINTINV[0]
print '1� Autovetor ='
print RESINTINV[1]
print 'N�mero de itera��es =', RESINTINV[2]
#Limpa a lista RESINTINV
RESINTINV=[]
# Inicia os autovalores
lbdkmais1, lbd = 0, 0
# Lista de coeficientes cj
listcj=[]
################## Ortonormaliza��o por Gram-Schmidt ##################
for u in range(ap-1):
    # Inicia a vari�vel l�gica (a) e inicia o contador de itera��o
    a, n = 0, 0
    while a!= 1:
        if n>=1:
            lbd = lbdkmais1
        # C�lculo do Yk
        for i in range(m):
            soma =0
            for j in range(m):
                soma = soma + M[i][j]*Xotkmais1[j][0]
            Yk[i][0] = soma
        # C�lculo de Xk+1
        Xkmais1 = function_Givens.sl(K,Yk)
        # C�lculo de cj
        listcj=[]
        for j in range(len(listegvectors)):
            cj =0
            for i in range(m):
                soma =0
                for t in range(m):
                    soma = soma + M[i][t]*Xkmais1[t][0]
                #print soma
                cj = cj + listegvectors[j][i][0]*soma
            listcj.append(cj)
        # C�lculo de Xotkmais1
        for j in range(len(listcj)):
            cj=listcj[j]
            for i in range(m):
                Xkmais1[i][0] = Xkmais1[i][0] - cj*listegvectors[j][i][0]
                Xotkmais1[i][0] = Xkmais1[i][0]
        # C�lculo do autovalor pelo quociente de Rayleigh
        smM, smK = 0, 0
        for i in range(m):
            somaM, somaK = 0, 0
            for t in range(m):
                somaM = somaM + M[i][t]*Xotkmais1[t][0]
                somaK = somaK + K[i][t]*Xotkmais1[t][0]
            #print somaM, somaK
            smM = smM + Xotkmais1[i][0]*somaM
            smK = smK + Xotkmais1[i][0]*somaK
        # C�lculo do autovalor
        lbdkmais1 = smK/smM
        # Normaliza o autovetor pela massa
        for i in range(m):
            Xotkmais1[i][0] = Xotkmais1[i][0]/(smM**0.5)
        # Verifica a converg�ncia
        if (((lbdkmais1 - lbd)**2)**0.5)/((lbdkmais1**2)**0.5) <= tol or n>50:
            a=1
        n=n+1
    ## Armazena resultados ##
    # Armazena autovetores
    listegvectors.append(Xotkmais1)
    # Imprime resultados de cada autopar
    print str(u+2)+'� Autovalor =',lbdkmais1
    print str(u+2)+'� Autovetor ='
    print Xotkmais1
    print 'N�mero de itera��es =', n
    # Reinicia Xotkmais1    
    Xotkmais1 = ones([m,1],float)
