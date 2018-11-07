! A fortran95 program for G95
! By WQY

program main
  !implicit none
  !integer re_i
  !write(*,*) "Hello World!"
  !re_i = system("pause")

    implicit none
    real :: x0(2)
    real, dimension(:), allocatable :: x

    x0(1) = 0
    x0(2) = 2

    x = newton_raph(x0, 1e-12, .true., 100)
    write(*,*) x(1), x(2)

    write(*,*) x(1), x(2), x(3)

    !x0(1) = 1
    !x0(2) = 2
    !x = newton_raph(x0, 1e-12, .true., 100)


contains
    function trans(b) result(t)
        integer :: n, i
        real, dimension(:) :: b
        real, dimension(:,:), allocatable :: t
        n = size(b)
        allocate(t(n, 1))
        do i = 1,n
            t(i, 1) = b(i)
        end do
    end function trans

    function f(p) result(res)
        !sistema de equações
        !eq1: x + y ^ 2 = 4
        !eq2: e ^ x + xy = 3
        real, dimension(:) :: p
        real, dimension(:), allocatable :: res
        integer :: n
        n = size(p)

        res(1) = p(1)+(p(2)**2)-4
        res(2) = exp(p(1))+p(1)*p(2)-3
    end function f

    function jac(x) result(J)
        real :: dx=1e-10
        real, dimension(:) :: x
        real, dimension(:,:), allocatable :: J
        real, dimension(:), allocatable :: xn, dy, dif
        integer :: n, i, k

        n = size(x)
        allocate(J(n, n), xn(n), dy(n), dif(n))
        do k = 1,n
            xn = x
            x(k) = xn(k)+dx
            dy = f(xn) - f(x)
            dif = dy/dx
            do i = 1,n
                J(i, k) = dif(i)
            end do
        end do
    end function jac

    function solve(M, b) result(x)
        real, dimension(:) :: b, M(:, :)
        real, dimension(:), allocatable :: x
        real :: lambda=1, es=1e-6, dummy, soma, old, ea
        integer :: n, iter, sentinela, i, j

        n = size(b)
        allocate(x(n))

        do i = 1,n
            x(i) = 0
        end do

        do i = 1,n
            dummy = M(i, i)
            do j = 1,n
                M(i, j) = M(i, j)/dummy
            end do
            b(i) = b(i)/dummy
        end do

        do i = 1,n
            soma = b(i)
            do j = 1,n
                if (i/=j) then
                    soma = soma - M(i, j)*x(j)
                end if
            end do
            x(i) = soma
        end do

        iter = 1
        sentinela = 1

        do while (sentinela==1)
            do i = 1,n
                old = x(i)
                soma = b(i)
                do j = 1,n
                    if (i/=j) then
                        soma = soma - M(i, j)*x(j)
                    end if
                end do
                x(i) = lambda*soma + (1-lambda)*old

                if ((sentinela==1).and.(x(i)/=0)) then
                    ea = abs((x(i)-old)/x(i))*100
                end if
                if (ea>es) then
                    sentinela = 0
                end if
            end do
        end do
    end function solve

    function newton_raph(x0, tol, iter, n_tot) result(x)
        real, dimension(:) :: x0
        real, dimension(:), allocatable :: x, e, F, S
        real :: tol
        integer :: n_tot, n0, n=0, i, j, k
        logical :: iter

        n0 = size(x0)
        tol = abs(tol)
        n_tot = abs(n_tot)
        allocate(x(n0), e(n0), F(n0), S(n0))
        do i = 1,n0
            x(i) = real(x0(i))
        end do
        e = maxval(abs(f(x)))

        do while ((e>tol).and.(n<=n_tot))
            n = n+1
            F = f(x)
            e = maxval(abs(F))
            J = jac(x)
            if (size(F)==1) then
                x = x - (F/J)
            else
                S = solve(J, -F)
                x = x+S
            end if
        end do

        if (iter .eqv. .true.) then
            write(*,*) "Total de iterações: ", n
        end if
        if (n>=n_tot) then
            write(*,*) "Processo parou, número de iterações limite atingido", n
        end if
    end function newton_raph
end program main
