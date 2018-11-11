function LU(matriz, vetor_b)
    n = length(vetor_b)  # Give us total of lines
    x = zeros(n)
    # (2) Fill L matrix and its diagonal with 1
    L = zeros(n, n)
    for i = 1:n
        L[i, i] = 1
    end
    # (3) Fill U matrix
    U=zeros(n,n)
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
    #print(L)
    y = zeros(n)
    # (5) Perform substitution Ly=b
    for i = 1:n
        y[i] = vetor_b[i] / L[i, i]
        
        for k = 1:i-1
            y[i] = y[i]  - y[k]*L[i, k]
            #print(y[i], "\n")
        end
        #print("y[", i, "] =", y[i], "\n")
    end
    #print(y, "\n")
    n = length(y)
    #print(n)
    x = copy(y)
    #print("\n", x, "\n")
    # (6) Perform substitution Ux=y
    for i = (n:-1:1)
    
        for k = 1+i:n
            x[i] = x[i] - x[k]*U[i, k]
        end
        x[i] = x[i]/U[i, i]
    end
    #print(x, "\n")

    return x
end

function f(p)
    # sistema de equações
    # eq1: x + y ^ 2 = 4
    # eq2: e ^ x + xy = 3
    a = p[2] + p[3] - exp(-p[1])
    b = p[1] + p[3] - exp(-p[3])
    c = p[1] + p[2] - exp(-p[3])
    return [a c b]
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


    x0 = [0.5, 1, 5]
    x = newton_raph(x0, 1e-12, true, 100)
    print("\nSolução= ",x, "\n")
