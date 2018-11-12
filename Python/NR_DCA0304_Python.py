import numpy as np
import time
from numpy import sin, exp
from math import pi

#DECLARACAO DAS FUNCOES
def func(p):
    '''
    Sistema de equacoes

        eq1: 1/2*sin(x1*x2)-(x2/(4*pi))-(x1/2),
        eq2: (1-1/(4*pi))*((exp(2*x1))-exp(1))-((exp(1)*x2)/pi)-2*exp(1)*x1)

        eq3: x2+x3-exp(-x1)
        eq4: x1+x3-exp(-x3)
        eq5: x1+x2-exp(-x3)
    '''

    x1, x2 = p

    return np.array((1/2*sin(x1*x2)-(x2/(4*pi))-(x1/2),
                     (1 - 1 / (4 * pi)) * ((exp(2 * x1)) - exp(1)) - ((exp(1) * x2) / pi) - 2 * exp(1) * x1))


#MATRIZ JACOBIANA
def jac(f, x, dx=1e-10):
    x = np.array(x)
    n = len(x)
    J = np.zeros((n, n))
    for j in range(n):          #Calculo numerico da matriz jacobiana
        xn = np.copy(x)
        xn[j] = xn[j] + dx
        dy = np.array(f(xn)) - np.array(f(x))
        dif = dy / dx
        for i in range(n):
            J[i, j] = dif[i]
    return J


# METODO DA FATORACAO LU
def LU(matriz, vetor_b):

    n = len(matriz)

    L = [[0 for i in range(n)] for i in range(n)]
    for i in range(n):
        L[i][i] = 1         #Preenche L com a Matriz Identidade

    U = [[0 for i in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(n):
            U[i][j] = matriz[i][j]       #Atribui a Matriz J a Matriz U


    # Encontra as Matrizes U e L
    for i in range(n-1):
        for k in range(i + 1, n):
            c = U[k][i] / U[i][i]
            L[k][i] = c  # Armazena o multiplicador
            for j in range(n):
                U[k][j] -= c * U[i][j]  # Multiplica com o pivo da linha e  subtrai


        for k in range(i + 1, n):
            U[k][i] = 0


    # Resolve o Sistema Ly=b
    y = [0 for i in range(n)]

    for i in range(0, n, 1):
        y[i] = vetor_b[i] / L[i][i]
        for k in range(0, i, 1):
            y[i] -= y[k] * L[i][k]


    # Resolve o Sistema Ux=y
    x = [y[i] for i in range(n)]

    for i in range(n-1, -1, -1):
        for k in range(i+1, n):
            x[i] -= x[k] * U[i][k]
        x[i] /= U[i][i]
    return x

#METODO NEWTON_RAPHSON
def newton_raph(func, x0, tol, iter):

    #DECLARACAO DAS VARIAVEIS
    tol = abs(tol)      #Modulos de tol e iter
    iter = abs(iter)
    x = (np.array(x0)).astype(np.float)    #Cria o vetor x
    e = max(abs(np.array(func(x))))  #Erro
    n = 0      #Atribue zero a variavel de interacoes totais

    #Print da interacao 0
    if iter == 0:
        print('\nTotal de Iteracoes: ' + str(n))
        return x

    else:

        while (e > tol and n <= iter):
            n += 1
            F = np.array(func(x))
            e = max(abs(F))
            J = jac(func, x)
            if len(F) == 1:
                x = x - F / J
            else:

                S = LU(J, -F)
                x = x + S

        print('\nTotal de Iteracoes: ' + str(n))
        if n >= iter:  #Condicao de parada
            print("PROCESSO PAROU, numero de iteracoes limite atingido!")
            return x
        else:
            return x


if __name__ == "__main__":

    inicio = time.time()
    x0 = [0.5, -2]  #Chute Inicial
    iter = 100      #Quantidade de Interacoes maximas
    tol = 1e-12     #tolerancia

    x = newton_raph(func, x0, tol, iter)   #Chama o metodo de Newton_Raphson

    print('Solucao: {} ' .format(x))

    fim = time.time()
    print('Tempo gasto: {}'.format(fim - inicio))
