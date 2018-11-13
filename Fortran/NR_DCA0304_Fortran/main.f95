program main
    implicit none
    double precision :: x0(3), x(size(x0)), tol, iter
    !Chute inicial
    x0(1) = 0.5
    x0(2) = 0.5
    x0(3) = 0.5
    iter = 100 !Numero de iteracoes
    tol = 1e-12 !Tolerancia

    x = newton_raph(x0, tol, .true., iter)
    write(*,*) x(1), x(2), x(3) !Solucao



contains
    
    !DECLARACAO DAS FUNCOES
    function f(p) result(res)
        IMPLICIT NONE
        double precision :: p(:), res(size(p))
        
        !Sistema de equacoes

        !eq1: 1/2*sin(x1*x2)-(x2/(4*pi))-(x1/2),
        !eq2: (1-1/(4*pi))*((exp(2*x1))-exp(1))-((exp(1)*x2)/pi)-2*exp(1)*x1)

        !eq3: x2+x3-exp(-x1)
        !eq4: x1+x3-exp(-x3)
        !eq5: x1+x2-exp(-x3)

        res(1) = p(2) + p(3) - exp(-p(1))
        res(2) = p(1) + p(3) - exp(-p(3))
        res(3) = p(1) + p(2) - exp(-p(3))
    end function f
    !MATRIZ JACOBIANA
    function jac(x) result(jc)
        IMPLICIT NONE
        double precision :: dx=1e-12, x(:), jc(size(x), size(x)), xn(size(x)), dy(size(x)), dif(size(x))
        integer :: n, i, k

        n = size(x)
        do k = 1,n !Calculo numerico da matriz jacobiana
            xn = x
            xn(k) = xn(k)+dx
            dy = f(xn) - f(x)
            dif = dy/dx
            do i = 1,n
                jc(i, k) = dif(i)
            end do
        end do
    end function jac
    ! METODO DA FATORACAO LU
    function LU(M, b) result(x)
        IMPLICIT NONE
        double precision :: b(:), M(:, :), x(size(b)), y(size(b)), U(size(b), size(b)), L(size(b), size(b)), c
        integer :: n, i, j, k

        n = size(b)

        do i = 1,n
            do j = 1,n
                L(i, j) = 0
                U(i, j) = 0
            end do
        end do
        do i = 1,n
            x(i) = 0
            y(i) = 0
            L(i, i) = 1 !Preenche L com a Matriz Identidade
        end do

        U = M !Atribui a Matriz J a Matriz U

        do i = 1,(n-1)
            do k = (i+1),n
                c = U(i, k) / U(i, i)
                L(k, i) = c !Armazena o multiplicador
                do j = 1,n
                    U(k, j) = U(k, j) - c*U(i, j) !Multiplica com o pivo da linha e subtrai
end
                end do
            end do
            do k = (i+1),n
                U(k, i) = 0
            end do
        end do
        
        !Resolve o Sistema Ly=b
        do i = 1,n
            y(i) = b(i) / L(i, i)
            do k = 1,(i-1)
                y(i) = y(i) - y(k)*L(i, k)
            end do
        end do

        x = y
        !Resolve o Sistema Ux=y
        do i = n,1,-1
            do k = (i+1),n
                x(i) = x(i) - x(k)*U(i, k)
            end do
            x(i) = x(i)/U(i, i)
        end do
    end function LU
    !METODO NEWTON_RAPHSON
    function newton_raph(x0, tol, iter, n_tot) result(x)
        IMPLICIT NONE
        double precision :: x0(:), x(size(x0)), tol, e, func(size(x0)), S(size(x0)), jc(size(x0), size(x0))
        integer :: n_tot, n0, n=0, i
        logical :: iter

        n0 = size(x0)

        do i = 1,n0
            x(i) = x0(i)
        end do
        e = maxval(abs(f(x)))

        do while ((e>tol).and.(n<=n_tot))
            n = n+1
            func = f(x)
            e = maxval(abs(func))
            jc = jac(x)
            if (size(func)==1) then
                x = x - (func(1)/jc(1, 1))
            else
                S = LU(jc, -func)
                x = x+S
            end if
        end do

        if (iter .eqv. .true.) then
            write(*,*) "Total de iteracoes: ", n
        end if
        if (n>=n_tot) then
            write(*,*) "Processo parou, numero de iteracoes limite atingido", n
        end if

    end function newton_raph

end program main
