import numpy as np
import time
from numpy import sin, exp
from math import pi

def f(p):
    '''
    Sistema de equacoes

        eq1: 1/2*sin(x*y)-(y/(4*pi))-(x/2),
        eq2: (1-1/(4*pi))*((exp(2*x))-exp(1))-((exp(1)*y)/pi)-2*exp(1)*x)

        eq3: y+z-exp(-x)
        eq4: x+y-exp(-z)
        eq5: x+z-exp(-z)
    '''

    x, y = p

    return np.array((1/2*sin(x*y)-(y/(4*pi))-(x/2),
                     (1 - 1 / (4 * pi)) * ((exp(2 * x)) - exp(1)) - ((exp(1) * y) / pi) - 2 * exp(1) * x))

def jac(f, x, dx=1e-10):
    x = np.array(x)
    n = len(x)
    J = np.zeros((n, n))
    for j in range(n):
        xn = np.copy(x)
        xn[j] = xn[j] + dx
        dy = np.array(f(xn)) - np.array(f(x))
        dif = dy / dx
        for i in range(n):
            J[i, j] = dif[i]
    return J

def LU(matriz, vetor_b):
    n = len(matriz)  

    L = [[0 for i in range(n)] for i in range(n)]
    for i in range(0, n):
        L[i][i] = 1

    U = [[0 for i in range(0, n)] for i in range(n)]
    for i in range(0, n):
        for j in range(0, n):
            U[i][j] = matriz[i][j]
    n = len(U)


    for i in range(0, n-1): 
        for k in range(i + 1, n):
            c = U[k][i] / U[i][i]
            L[k][i] = c  
            for j in range(0, n):
                U[k][j] -= c * U[i][j] 

        for k in range(i + 1, n):
            U[k][i] = 0
    n = len(L)

    y = [0 for i in range(n)]
    for i in range(0, n, 1):
        y[i] = vetor_b[i] / L[i][i]
        for k in range(0, i, 1):
            y[i] -= y[k] * L[i][k]

    n = len(U)

    x = [y[i] for i in range(0, n)]
    for i in range(n-1, -1, -1):
        for k in range(i+1, n):
            x[i] -= x[k] * U[i][k]
        x[i] /= U[i][i]
    return x






def newton_raph(func, x0, tol, iter):
    """
    -----------------------------------------------------------------------------------------------
	    Parametros:

	    func:       //E  uma funcao que deve conter todas as equacoes do sistema.
	    x0:         //Uma lista ou array com chute inicial para cada variavel.
	    tol:        //Determina a tolerancia minima dos valores da matriz F como criterio de parada.
	                    Este valor deve ser proximo de zero, logo e aceitavel tol = 1e-12.
	    iter:      //Representa um numero de iteracoes maximo como criterio de parada.
    -----------------------------------------------------------------------------------------------
    """
    tol = abs(tol)
    iter = abs(iter)
    x = (np.array(x0)).astype(np.float)
    e = max(abs(np.array(func(x))))
    n = 0
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
        if n >= iter:
            print("PROCESSO PAROU, numero de iteracoes limite atingido!")
            return x
        else:
            return x


if __name__ == "__main__":

    x0 = [0.5, 0.5]
    iter = 100
    tol = 1e-12

    x = newton_raph(f, x0, tol, iter)

    print('Solucao: {} ' .format(x))

    fim_sis = time.process_time()
    print('Tempo gasto: {}'.format(fim_sis))
