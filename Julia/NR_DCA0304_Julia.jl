

function solve_sis(M, b)
    n = length(b)
    x = zeros(n)
    xo = zeros(n)
    del=1e-6
    intera = 10
    k = 0
    test = 1

    while (k<intera)&(test>del)
        k = k+1
        xo = x #Jacobi
        for i = 1:n
            summ = 0;
            for j = 1:n
                summ = summ + M[i, j]*xo[j]*(i!=j)
            end
            x[i] = (b[i]-summ)/M[i, i]
            test = x[i]-xo[i]
        end
    end
    return x
end

function LU(matriz, vetor_b)
    n = length(vetor_b)  # Give us total of lines
    y = zeros(n)
    x = zeros(n)
    # (2) Fill L matrix and its diagonal with 1
    L = zeros(n, n)
    for i = 1:n
        L[i, i] = 1
    end

    # (3) Fill U matrix
    U = matriz
    
    
    # (4) Find both U and L matrices
    for i = 1:n-1
        for k = i+1:n
            c = U[k, i] / U[i, i]
            L[k, i] = c # (4.4) Store the multiplier
            for j = 1:n
                U[k, j] = U[k, j] - c*U[i, j] # Multiply with the pivot line and subtract
            end
        end
        # (4.5) Make the rows bellow this one zero in the current column
        for k = i+1:n
            U[k, i] = 0
        end
    end
        
    # (5) Perform substitution Ly=b
    for i = 1:n
        y[i] = vetor_b[i] / L[i, i]
        for k = 1:i
            y[i] = y[i] - y[k]*L[i, k]
        end
        #print("y[", i, "] =", y[i], "\n")
    end
    x = copy(y)
    #print("\n", x, "\n")
    # (6) Perform substitution Ux=y
    for i = n-1:-1:1
        for k = i+1:n
            x[i] = x[i] - x[k]*U[i, k]
        end
        x[i] = x[i]/U[i, i]
    end


    return x
end

function f(p)
    # sistema de equações
    # eq1: x + y ^ 2 = 4
    # eq2: e ^ x + xy = 3
    a = p[2] + p[3] - exp(-p[1])
    b = p[1] + p[3] - exp(-p[3])
    c = p[1] + p[2] - exp(-p[3])
    return [a b c]
end

function jac(x, dx=1e-10)
    n = length(x)
    J = zeros(n, n)
    for j = 1:n
        xn = copy(x)
        xn[j] = xn[j] +dx
        dy = f(xn) - f(x)
        dif = dy/dx
        for i = 1:n
            J[i, j] = dif[i]
        end
    end
    return J
end

function newton_raph(x0, tol, iter, n_tot)
    tol = abs.(tol)
    n_tot = abs.(n_tot)
    x = convert(Array{Float64}, x0)
    e = maximum(abs.(f(x)))
    n = 0
    while (e>tol)&(n<=n_tot)
        n += 1
        F = f(x)
        e = maximum(abs.(F))
        J = jac(x)
        if length(F) == 1
            x = x - F / J
        else
            S = LU(J, -F)
            x = x+S
        end
    end
    if iter==true
        print("Total de Iterações: ", string(n))
    end
    if n>=n_tot
        print("Processo parou, número de iterações limite atingido")
    else
        return x
    end
end


    x0 = [0.5, 0.5, 0.5]
    x = newton_raph(x0, 1e-6, true, 1)
    print(x)
