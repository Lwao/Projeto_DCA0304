function solve(M, b, prec)
    n = length(b)
    x0 = zeros(n)
    x0 = convert(Array{Float64}, x0)
    x = zeros(n)
    x = convert(Array{Float64}, x)
    intera = 50

    for k = 1:intera
        x0 = x #jacobi
        for i = 1:n
            summ = float(0)
            for j = 1:n
                summ = summ + M[i, j]*x0[j]*(i!=j)
            end
            x[i] = (1/M[i, i])*(b[i]-summ)
            #x0 = x #gauus-siedel
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
            S = solve(J, -F, tol)
            x = x + S
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


    x0 = [0, 2]
    x = newton_raph(x0, 1e-12, true, 100)
    print("\nSolução: {", string(x[1]), "}")
