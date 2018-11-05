import numpy as np


def f(p):
    '''
    sistema de equações
        eq1: x + y ^ 2 = 4
        eq2: e ^ x + xy = 3
    '''
    x, y = p
    return np.array((x + (y ** 2) - 4,
                     np.exp(x) + x * y - 3))

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


def newton_raph(func, x0, tol, iter, n_tot):
    """
    -----------------------------------------------------------------------------------------------
	    Esta função implementa o método Newton-Raphson para sistemas de equações não lineares.
	    Parâmetros:

	    func:       //É  uma função que deve conter todas as equações do sistema.
	    x0:         //Uma lista ou array com chute inicial para cada variável.
	    tol:        //Determina a tolerância mínima dos valores da matriz F como critério de parada.
	             Este valor deve ser próximo de zero, logo é aceitável tol = 1e-12.
	    iter:       //Uma chave para que seja exibido o total de iterações no fim do processo.
	    n_tot:      //Representa um número de iterações máximo como critério de parada.
    -----------------------------------------------------------------------------------------------
    """
    tol = abs(tol)
    n_tot = abs(n_tot)
    x = (np.array(x0)).astype(np.float)
    e = max(abs(np.array(func(x))))
    n = 0
    while (e > tol and n <= n_tot):
        n += 1
        F = np.array(func(x))
        e = max(abs(F))
        J = jac(func, x)
        if len(F) == 1:
            x = x - F / J
        else:
            S = np.linalg.solve(J, -F)
            x = x + S
    if iter == True:
        print('Total de Iterações: ' + str(n))
    if n >= n_tot:
        print("Processo parou, número de iterações limite atingido")
    else:
        return x

if __name__ == "__main__":

    x0 = [0, 2]
    x = newton_raph(f, x0, 1e-12, True, 100)
    print('Solução: {} ' .format(x))