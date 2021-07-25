# -*- coding: cp1252 -*-

################################################################################
#            Integra��o num�rica pela regra trapezoidal composta               #
################################################################################
from numpy import array, zeros
from copy import *
import function_NDDIP

#################################### Entrada ###################################
# Lista de coordenadas x (listx) e F(x) (listFx) para interpola��o
listx=array([-250,-200,-100,0,100,300])
listFx=array([0.0163,0.318,0.699,0.87,0.941,1.04])
# Intervalo de integra��o
a, b = 0, 300
# N�mero de pontos de integra��o
m=3
################################################################################
print 'N�mero de pontos de integra��o ='
print m
### Define os pontos de integra��o ###
# Vetores de coordenadas x e Fx
vpix, vpifx = zeros([m,1],float), zeros([m,1],float)
vpix[0][0], vpix[m-1][0] = a, b
Pt = a
# Define intervalos entre os pontos igualmente espa�ados
for i in range(1,m-1):
    Pt = Pt + (b-a)/(m-1)*1.0
    vpix[i][0] = Pt
# Coordenadas Fx
vpifx = function_NDDIP.NI(listx,listFx,vpix)
print 'Coordenadas x dos pontos de integra��o ='
print vpix
print 'Coordenadas Fx dos pontos de integra��o ='
print vpifx
# Inicia o valor da integral no intervalo
H=0
# C�lculo da integral no intervalo
for i in range(1,m):
    h = (vpix[i][0]-vpix[i-1][0])
    H = H + (h/2)*(vpifx[i][0]+vpifx[i-1][0])
print 'Valor num�rico da integral H =', H
