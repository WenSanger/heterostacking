subroutine twist_angle(m1,n1,m2,n2,angle)
    implicit none
    real(kind=8)::angle, cos_theta
    integer :: m1,n1,m2,n2
    real(kind=8),parameter::pi = 3.141592653589793238
    cos_theta = (m1*m2 + n1*n2 + 0.5d0*m1*n2 + 0.5d0*m2*n1) / (sqrt((m1**2+n1**2+m1*n1)*1.0d0)*sqrt((m2**2+n2**2+m2*n2)*1.0d0))
    angle = acos(cos_theta)/pi*180
end subroutine
