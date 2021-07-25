# -*- coding: cp1252 -*-

################################################################################
#              Integração numérica pela regra Simpson composta                 #
################################################################################
from numpy import array, zeros
from copy import *
import function_NDDIP

#################################### Entrada ###################################
# Lista de coordenadas x (listx) e F(x) (listFx) para interpolação
listx=array([-250,-200,-100,0,100,300])
listFx=array([0.0163,0.318,0.699,0.87,0.941,1.04])
# Intervalo de integração
a, b = 0, 300
# Número de pontos de integração
m = 3 # sempre números impares (3, 5, 7,...)
################################################################################
print 'Número de pontos de integração ='
print m
### Define os pontos de integração ###
# Vetores de coordenadas x e Fx
vpix, vpifx = zeros([m,1],float), zeros([m,1],float)
vpix[0][0], vpix[m-1][0] = a, b
Pt = a
# Define intervalos entre os pontos igualmente espaçados
h = (b-a)/(m-1)*1.0
print 'Valor de h =',h
for i in range(1,m-1):
    Pt = Pt + h
    vpix[i][0] = Pt
# Coordenadas Fx
vpifx = function_NDDIP.NI(listx,listFx,vpix)
print 'Coordenadas x dos pontos de integração ='
print vpix
print 'Coordenadas Fx dos pontos de integração ='
print vpifx
# Inicia o valor da integral no intervalo
H=0
# Cálculo da integral de Simpson no intervalo
H = vpifx[0][0]+vpifx[m-1][0]
for i in range(1,m-1):
    # verifica se i é par ou ímpar
    if i%2==0:
        H = H + (2*vpifx[i][0])
    else:
        H = H + (4*vpifx[i][0])
H = (h/3)*H
print 'Valor numérico da integral H =', H
