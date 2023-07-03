program main
!program calculating parameter of the periodic structure in MoSe2/WSe2 twisted heterostructure. 
    implicit none
    real(kind=8) :: a1, a2, theta, f, f_limit, percent, biggest
    integer :: m1, n1, m2, n2 
    f_limit = 0.01d0!A, upper limit of distance_func
    a2 = 3.315d0!3.2890000343
    a1 = 3.181d0!3.2820000648 
    open(11,file = "./result.dat")
    write(11,"(A,4X,A,4X,A,4X,A,4X,A,4X,A)")"theta","m1","n1","m2","n2","point_distance"
    biggest = 0.0d0
    !$omp parallel do
    do m1 = 0, 50!300
        do n1 = 0, m1
            do m2 = 0, 50!300
                do n2 = 0, m2
                    call distance_func(m1, n1, a1, m2, n2, a2, f, percent)
                    if (f .lt. f_limit) then
                        call twist_angle(m1,n1,m2,n2,theta)
                        !write(11,"(F13.10,4X,I3,4X,I3,4X,I3,4X,I3,4X,F12.10,4x,F12.10)") theta, m1, n1, m2, n2, f, percent
                        if (theta .ge. 0.064 ) then
                            write(11,"(F13.10,4X,I3,4X,I3,4X,I3,4X,I3,4X,F12.10)") theta, m1, n1, m2, n2, f
                        end if
                        if (percent .gt. biggest) then
                            biggest = percent
                        end if
                    end if
                end do
            end do
        end do
    end do
    !$omp end parallel do
    write(11, *) "biggest percentage = ", biggest
end program main
