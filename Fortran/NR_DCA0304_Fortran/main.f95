! A fortran95 program for G95
! By WQY

program main
    implicit none
    double precision :: x0(3), x(size(x0)), tol=1e-12

    x0(1) = 0.5
    x0(2) = 0.5
    x0(3) = 0.5

    x = newton_raph(x0, tol, .true., 100)
    write(*,*) x(1), x(2), x(3)



contains

    function f(p) result(res)
        IMPLICIT NONE
        double precision :: p(:), res(size(p))

        res(1) = p(2) + p(3) - exp(-p(1))
        res(2) = p(1) + p(3) - exp(-p(3))
        res(3) = p(1) + p(2) - exp(-p(3))
    end function f

    function jac(x) result(jc)
        IMPLICIT NONE
        double precision :: dx=1e-12, x(:), jc(size(x), size(x)), xn(size(x)), dy(size(x)), dif(size(x))
        integer :: n, i, k

        n = size(x)
        do k = 1,n
            xn = x
            xn(k) = xn(k)+dx
            dy = f(xn) - f(x)
            dif = dy/dx
            do i = 1,n
                jc(i, k) = dif(i)
            end do
        end do
    end function jac

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
            L(i, i) = 1
        end do

        U = M

        do i = 1,(n-1)
            do k = (i+1),n
                c = U(i, k) / U(i, i)
                L(k, i) = c
                do j = 1,n
                    U(k, j) = U(k, j) - c*U(i, j)
                end do
            end do
            do k = (i+1),n
                U(k, i) = 0
            end do
        end do

        do i = 1,n
            y(i) = b(i) / L(i, i)
            do k = 1,(i-1)
                y(i) = y(i) - y(k)*L(i, k)
            end do
        end do

        x = y

        do i = n,1,-1
            do k = (i+1),n
                x(i) = x(i) - x(k)*U(i, k)
            end do
            x(i) = x(i)/U(i, i)
        end do
    end function LU

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
            write(*,*) "Total de iterações: ", n
        end if
        if (n>=n_tot) then
            write(*,*) "Processo parou, número de iterações limite atingido", n
        end if

    end function newton_raph

end program main
