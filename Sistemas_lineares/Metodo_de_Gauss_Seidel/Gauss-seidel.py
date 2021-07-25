# -*- coding: cp1252 -*-
# Resolu��o de um sistema linear ([K]{U}={R}) pelo m�todo de Gauss-Seidel
print 'Resolu��o de um sistema linear ([K]{U}={R}) pelo m�todo de de Gauss-Seidel'
from numpy import *

# Matriz [A]
K = array([[4,-1,-1,0,-1,0],[-1,4,-1,-1,0,0],[-1,-1,4,-1,-1,-1],[0,-1,-1,4,-1,0],[-1,0,-1,-1,4,-1],[0,0,-1,0,-1,4]],float)
#K = array([[5,-4,1,0],[-4,6,-4,1],[1,-4,6,-4],[0,1,-4,5]],float)
print 'Matriz [K] = '
print K

# Vetor {R}
R=array([4,3,3,3,3,1],float)
#R=array([0,1,0,0],float)
print 'Vetor {R} = '
print R

# Fator de relaxa��o
Beta=1.4

# Ordem da matriz
n=len(K) # conta o n�mero de colunas da matriz [K]

# Inicia vetores
Us=zeros([n,1],float)
DeltaU=zeros([n,1],float)
Us1=zeros([n,1],float) # vetor inicial
print 'Uo = '
print Us1


# Inicia vari�veis
Epsilon, SS, SS2, iteracao = 1, 0, 0, 0

while Epsilon > 10**(-5):
    for k in range(n):
        # Obt�m DeltaU correspondente ao DeltaU={R}-[K](L)*{U}(s+1), onde [K](L)
        #� a matriz triangular inferior e {U}(s+1) � o vetor da pr�xima itera��o.
        for l in range(1,n):
            SS=0
            for i in range(l):
                SS=SS + (K[l][i]*Us1[i])
            DeltaU[l]=R[l]-SS
        DeltaU[0]=R[0]
        # Obt�m DeltaU correspondente ao DeltaU = DeltaU-[K]T(L)*{U}(s), onde [K]T(L)
        #� a matriz triangular inferior transposta e {U}(s) � o vetor da itera��o corrente.
        for i in range(n-1):
            SS=0
            for l in range(i+1,n):
                DeltaU[i]=DeltaU[i]-(K[i][l]*Us[l])

        # Obt�m DeltaU correspondente ao DeltaU = beta*[K]-1(D)*(DeltaU-[K](D)*{U}(s)),
        #onde [K](D) � a matriz diagonal e {U}(s) � o vetor da itera��o corrente.
        for i in range(n):
            DeltaU[i]=DeltaU[i]-(K[i][i]*Us[i])
            DeltaU[i]=Beta*((K[i][i])**-1)*DeltaU[i]
            # Obt�m o Us1 da pr�xima itera��o
            Us1[i] = Us[i] + DeltaU[i]


    # Conta o n�mero de itera��es
    iteracao = iteracao + 1

    # Verifica a toler�ncia de converg�ncia
    SS, SS2 = 0, 0
    for j in range(n):
        SS = SS + (Us1[j] - Us[j])**2
        SS2 = SS2 + (Us1[j])**2
        
    Epsilon = (SS**0.5) / (SS2**0.5)
    print Us
    # Iguala o vetor {U}(s+1) ao {U}(s)
    for j in range(n):
        Us[j] = Us1[j]

print Us1
print iteracao











    

                
            
