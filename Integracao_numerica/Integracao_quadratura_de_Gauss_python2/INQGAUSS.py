# -*- coding: cp1252 -*-
################################################################################
#                   Integra��o num�rica pela quadratura de gauss               #
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
# Lista de coeficientes Cn para integra��o de tr�s pontos
listC=array([(5./9.),(8./9.),(5./9.)])
# Coordenadas dos tr�s pontos de integra��o de gauss
vpix=array([[-(0.6)**0.5],[0],[(0.6)**0.5]])
################################################################################
print 'N�mero de pontos de integra��o ='
m=len(vpix)
print m
### Define os pontos de integra��o ###
# Vetores de coordenadas x e Fx
vpifx = zeros([m,1],float)

# Transforma��o de coordenadas
A, B = (b-a)/2, (b+a)/2
# Obten��o da coordenadas x e Fx
for i in range(m):
    vpix[i][0] = (A*vpix[i][0]+B)
print 'Coordenadas x dos pontos de integra��o ='
print vpix
# Coordenadas Fx
vpifx = function_NDDIP.NI(listx,listFx,vpix)
#Inicia o valor da integral de Gauss
H=0.0
for i in range(m):
    vpifx[i][0] = (A*vpifx[i][0])
    H = H + vpifx[i][0]*listC[i]
print 'Coordenadas Fx dos pontos de integra��o ='
print vpifx
print 'Valor num�rico da integral H =', H
