import numpy as np
import time
from numpy import sin, exp
from math import pi


def f(p):
    '''
    Sistema de equações

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
    #print(matriz)
    #print(vetor_b)
    n = len(matriz)  # Give us total of lines

    # (2) Fill L matrix and its diagonal with 1
    L = [[0 for i in range(n)] for i in range(n)]
    for i in range(0, n):
        L[i][i] = 1
    #print(L)

    # (3) Fill U matrix
    U = [[0 for i in range(0, n)] for i in range(n)]
    for i in range(0, n):
        for j in range(0, n):
            U[i][j] = matriz[i][j]
    #print(U)
    n = len(U)

    # (4) Find both U and L matrices
    for i in range(0, n-1):  # for i in [0,1,2,..,n]
        # (4.3) Subtract lines
        for k in range(i + 1, n):
            c = U[k][i] / U[i][i]
            L[k][i] = c  # (4.4) Store the multiplier
            for j in range(0, n):
                U[k][j] -= c * U[i][j]  # Multiply with the pivot line and subtract

        # (4.5) Make the rows bellow this one zero in the current column
        for k in range(i + 1, n):
            U[k][i] = 0
    #print('L = {}'.format(L))
    #print('U = {}'.format(U))
    #print(U)
    n = len(L)
    #print(vetor_b)
    # (5) Perform substitutioan Ly=b
    y = [0 for i in range(n)]
    for i in range(0, n, 1):
        #print(vetor_b[i])
        #print(L[i][i])
        y[i] = vetor_b[i] / L[i][i]
        for k in range(0, i, 1):
            y[i] -= y[k] * L[i][k]
            #print(y[i])

        #print('y[{}]={}'.format(i, y[i]))
    #print('y = {}'.format(y))
    n = len(U)

    # (6) Perform substitution Ux=y
    x = [y[i] for i in range(0, n)]
    for i in range(n-1, -1, -1):
        for k in range(i+1, n):
            x[i] -= x[k] * U[i][k]
        x[i] /= U[i][i]
    #print('x = {}'.format(x))
    return x



'''def gauss_seidel(matriz, vetor_b, iteracoes):

    dimensao = len(matriz)
    chute_inicial = dimensao*[0]

    x_anterior = [0.0 for i in range(dimensao)]
    for i in range(iteracoes):
        for j in range(dimensao):
            x_anterior[j] = chute_inicial[j]
        for j in range(dimensao):
            soma = 0.0
            for k in range(dimensao):
                if k != j:
                    soma += matriz[j][k] * chute_inicial[k]
            chute_inicial[j] = (vetor_b[j] - soma) / matriz[j][j]

    return chute_inicial'''



def newton_raph(func, x0, tol, iter):
    """
    -----------------------------------------------------------------------------------------------
	    Parâmetros:

	    func:       //É  uma função que deve conter todas as equações do sistema.
	    x0:         //Uma lista ou array com chute inicial para cada variável.
	    tol:        //Determina a tolerância mínima dos valores da matriz F como critério de parada.
	                    Este valor deve ser próximo de zero, logo é aceitável tol = 1e-12.
	    iter:      //Representa um número de iterações máximo como critério de parada.
    -----------------------------------------------------------------------------------------------
    """
    tol = abs(tol)
    iter = abs(iter)
    x = (np.array(x0)).astype(np.float)
    e = max(abs(np.array(func(x))))
    n = 0
    if iter == 0:
        print('\nTotal de Iterações: ' + str(n))
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
                #S = gauss_seidel(J, -F, iter)
                #S = np.linalg.solve(J, -F)
                x = x + S

        print('\nTotal de Iterações: ' + str(n))
        if n >= iter:
            print("PROCESSO PAROU, número de iterações limite atingido!")
            return x
        else:
            return x


if __name__ == "__main__":

    x0 = [0.5, 0.5]
    iter = 100
    tol = 1e-12

    x = newton_raph(f, x0, tol, iter)

    print('Solução: {} ' .format(x))

    fim_sis = time.process_time()
    print('Tempo gasto: {}'.format(fim_sis))
