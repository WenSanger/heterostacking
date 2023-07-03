subroutine distance_func(m1,n1,a1,m2,n2,a2,pnt_diff, percent)
    implicit none
    integer :: m1,n1,m2,n2
    real(kind=8) :: a1, a2, pnt_diff, percent
    pnt_diff = abs(a2*sqrt(m2**2.0d0+n2**2.0d0+m2*n2*1.0d0)-a1*sqrt(m1**2.0d0+n1**2.0d0+m1*n1*1.0d0))
    percent = pnt_diff / abs(a2*sqrt(m2**2.0d0+n2**2.0d0+m2*n2*1.0d0))
end subroutine
