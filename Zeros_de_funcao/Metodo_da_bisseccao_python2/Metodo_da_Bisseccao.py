# -*- coding: cp1252 -*-
############################### Método da Bissecção ############################
############ Encontra o zero da função f(x) pelo método da bissecção ###########

# Importando módulos
from math import log, ceil

# Insira a função:
def f(x):
    return x**3 - 10*x**2 + 5

# Insira os limites do intervalo a ser avaliado:
a, b = 0.6, 0.8

# Insira o valor da tolerância
epsilon = 1.0e-9

# Inicia as variáveis
x1, x2 = a, b
f1 = f(x1); f2 = f(x2)

# Verifica se o intervalo dado apresenta pelo menos uma raiz
if f1*f2 > 0.0:
    print 'Intervalo não é válido, pois não contém raizes'
else:
    # A função ceil arredonda para o maior e mais próximo número inteiro
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
    
