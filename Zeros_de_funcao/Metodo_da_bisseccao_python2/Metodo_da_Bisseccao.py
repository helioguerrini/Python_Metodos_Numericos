# -*- coding: cp1252 -*-
############################### M�todo da Bissec��o ############################
############ Encontra o zero da fun��o f(x) pelo m�todo da bissec��o ###########

# Importando m�dulos
from math import log, ceil

# Insira a fun��o:
def f(x):
    return x**3 - 10*x**2 + 5

# Insira os limites do intervalo a ser avaliado:
a, b = 0.6, 0.8

# Insira o valor da toler�ncia
epsilon = 1.0e-9

# Inicia as vari�veis
x1, x2 = a, b
f1 = f(x1); f2 = f(x2)

# Verifica se o intervalo dado apresenta pelo menos uma raiz
if f1*f2 > 0.0:
    print 'Intervalo n�o � v�lido, pois n�o cont�m raizes'
else:
    # A fun��o ceil arredonda para o maior e mais pr�ximo n�mero inteiro
    n = int(ceil(log(abs(x2-x1)/epsilon)/log(2.0)))
    for i in range(n):
        x3 = 0.5*(x1 + x2); f3 = f(x3)
        if f3 == 0.0:
            print x3
            break
        if f2*f3 < 0.0:
            x1 = x3; f1 = f3
        else:
            x2 = x3; f2 = f3

    print 'A raiz no intervalo vale x =', (x1 + x2)/2.0
    
