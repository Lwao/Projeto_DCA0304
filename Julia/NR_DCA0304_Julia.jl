# METODO DA FATORACAO LU
function LU(matriz, vetor_b)
    n = length(vetor_b)  
    x = zeros(n)
    L = zeros(n, n)
    for i = 1:n  #Preenche L com a Matriz Identidade
        L[i, i] = 1
    end
    U=zeros(n,n)
    U = matriz #Atribui a Matriz J a Matriz U

    for i = 1:n-1
        for k = i+1:n
            c = U[k, i] / U[i, i]
            L[k, i] = c #Armazena o multiplicador
            for j = 1:n
                U[k, j] = U[k, j] - c*U[i, j] # Multiplica com o pivo da linha e subtrai
            end
        end
        for k = i+1:n
            U[k, i] = 0
        end
    end
    
    # Resolve o Sistema Ly=b
    y = zeros(n)
    for i = 1:n
        y[i] = vetor_b[i] / L[i, i]
        
        for k = 1:i-1
            y[i] = y[i]  - y[k]*L[i, k]
        end
    end
    n = length(y)
    
    # Resolve o Sistema Ux=y
    x = copy(y)
    for i = (n:-1:1)
        for k = 1+i:n
            x[i] = x[i] - x[k]*U[i, k]
        end
        x[i] = x[i]/U[i, i]
    end

    return x
end


#DECLARACAO DAS FUNCOES
function f(p)
    #Sistema de equacoes

    #eq1: 1/2*sin(x1*x2)-(x2/(4*pi))-(x1/2),
    #eq2: (1-1/(4*pi))*((exp(2*x1))-exp(1))-((exp(1)*x2)/pi)-2*exp(1)*x1)

    #eq3: x2+x3-exp(-x1)
    #eq4: x1+x3-exp(-x3)
    #eq5: x1+x2-exp(-x3)
    a = p[2] + p[3] - exp(-p[1])
    b = p[1] + p[3] - exp(-p[3])
    c = p[1] + p[2] - exp(-p[3])
    return [a c b]
end
#MATRIZ JACOBIANA
function jac(x, dx=1e-10)
    n = length(x)
    J = zeros(n, n)
    for j = 1:n #Calculo numerico da matriz jacobiana
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
#METODO NEWTON_RAPHSON
function newton_raph(x0, tol, iter, n_tot)
    #DECLARACAO DAS VARIAVEIS
    tol = abs.(tol) #Modulos de tol e iter
    n_tot = abs.(n_tot)
    x = convert(Array{Float64}, x0) #Cria o vetor x
    e = maximum(abs.(f(x))) #Erro
    n = 0 #Atribue zero a variavel de interacoes totais

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
        print("Total de Iteracoes: ", string(n))
    end
    if n>=n_tot  #Condicao de parada
        print("Processo parou, numero de iteracoes limite atingido")
    else
        return x
    end
end


    x0 = [0.5, 1, 5] # Chute inicial
    iter  = 100 # Numero de iteracoes
    tol = 1e-12 # Tolerancia
    @timev x = newton_raph(x0, tol, true, iter)
    print("\nSolucao= ",x, "\n")
