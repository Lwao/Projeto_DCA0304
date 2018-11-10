function solve_sis(M, b)
    n = length(b)
    x = zeros(n)
    xo = zeros(n)
    del=1e-6
    intera = 50
    k = 0
    test = 0
    
    while (k<intera)|(abs.(test)>del)
        k = k+1
        #xo = x #Jacobi
        for i = 1:n
            summ = 0;
            for j = 1:n
                summ = summ + M[i, j]*xo[j]*(i!=j)
            end
            x[i] = (b[i]-summ)/M[i, i]
            test = x[i]-xo[i]
            xo = x #Gauss-Sidel
        end
    end
    return x
end

function f(p)
    # sistema de equações
    # eq1: x + y ^ 2 = 4
    # eq2: e ^ x + xy = 3
    a = p[1] + (p[2] ^ 2) - 4
    b = exp(p[1]) + p[1] * p[2] - 3
    return [a b]
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
            print("\n", J, "\n")
            print("\n", F, "\n")
            S = solve_sis(J, -F)
            print("\n", S, "\n")
            x =  x+S
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


    x0 = [2, 3]
    x = newton_raph(x0, 1e-6, true, 100)
    #print(x)
