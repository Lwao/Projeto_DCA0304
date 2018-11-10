! A fortran95 program for G95
! By WQY

program main
    implicit none
    real :: x0(3), x(size(x0))

    x0(1) = 0
    x0(2) = 1
    x0(3) = 2
    
    x = newton_raph(x0, 1e-12, .true., 100)
    write(*,*) x(1), x(2), x(3)



contains

    function f(p) result(res)
        IMPLICIT NONE
        real :: p(:), res(size(p))

        res(1) = p(2) + p(3) - exp(-p(1))
        res(2) = p(1) + p(2) - exp(-p(3))
        res(3) = p(1) + p(2) - exp(-p(3))
    end function f

    function jac(x) result(J)
        IMPLICIT NONE
        real :: dx=1e-10, x(:), J(size(x), size(x))
        real, dimension(:), allocatable :: xn, dy, dif
        integer :: n, i, k

        n = size(x)
        allocate(xn(n), dy(n), dif(n))
        do k = 1,n
            xn = x
            x(k) = xn(k)+dx
            dy = f(xn) - f(x)
            dif = dy/dx
            do i = 1,n
                J(i, k) = dif(i)
            end do
        end do
        
        deallocate(xn, dy, dif)
    end function jac

    function solve(M, b) result(x)
        IMPLICIT NONE
        real :: b(:), M(:, :), x(size(b)), xo(size(b)), test=1, del=1e-6, summ
        integer :: n, itera=10, i, j, k=0

        n = size(b)
        
        do while ((k<itera).and.(test>del))
            k = k+1
            do i = 1,n
                summ = 0;
                do j = 1,n
                    if (i/=j) then
                        summ = summ + M(i, j)*xo(j)
                    end if
                end do
                x(i) = (b(i)-summ)/M(i, i)
            end do
        end do
    end function solve

    function newton_raph(x0, tol, iter, n_tot) result(x)
        IMPLICIT NONE
        real :: x0(:), x(size(x0)), tol, e, temp(2)
        real, dimension(:), allocatable :: func, S, jc(:, :)
        integer :: n_tot, n0, n=0, i, j, k
        logical :: iter

        n0 = size(x0)
        tol = abs(tol)
        n_tot = abs(n_tot)
        allocate(jc(n0, n0), func(n0), S(n0))
        do i = 1,n0
            x(i) = real(x0(i))
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
                S = solve(jc, -func)
                x = x+S
            end if
        end do

        if (iter .eqv. .true.) then
            write(*,*) "Total de iterações: ", n
        end if
        if (n>=n_tot) then
            write(*,*) "Processo parou, número de iterações limite atingido", n
        end if
        
        deallocate(func, S, jc)
        
    end function newton_raph

end program main
