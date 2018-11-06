function trans(b)
    n = length(b)
    t = zeros(n, 1)
    for i = 1:n
        t[i, 1] = b[i]
    end
    return t
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
            S = inv(J)*trans(-F)
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
    print(x)
