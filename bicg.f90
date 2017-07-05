module bicg
  use constants
  use matrixmain
  implicit none

contains
  subroutine bicgsolve(matrixin,rhs,soln)
    class(matrix),intent(inout) :: matrixin !no longer class-bound subroutine
    real(kind=dp),dimension(matrixin%matrixsize) :: r_0,rhs,soln,r_old,r_new,r_hat,v_new,v_old,x_old,x_new,p_old,p_new, h_old, h_new,&
    & ss,tt
    real(kind=dp) :: err !can't give intent to variables not in argumentlist
    real(kind=dp) :: alpha,rho_old,rho_new, beta, omega_old,omega_new
    integer :: count
    v_old=0.0_dp
    p_old=0.0_dp

    call matrixin%multiply(soln,r_0)
    r_0=rhs-r_0
    r_old=r_0
    r_hat=r_0
    rho_old=1.0_dp
    omega_old=1.0_dp
    alpha=1.0_dp
    x_old=soln
do count=1,10000
  print*, count, "===================="
!1

    rho_new=dot_product(r_hat,r_old)
!2
    beta=(rho_new/rho_old)*(alpha/omega_old)
!3

    p_new=r_old+beta*(p_old-omega_old*v_old)
!4
    call matrixin%multiply(p_new,v_new)

!5
    alpha=rho_new/(dot_product(r_hat,v_new))
!6

    h_new=x_old+alpha*p_new
!7
    if (all(abs(h_new-h_old)<=h_new*1.0E-10_dp)) then !h is accurate enough
      soln=h_new
      return
    end if
!8
      ss=r_old-alpha*v_new
!9
      call matrixin%multiply(ss,tt)
!10
      omega_new=(dot_product(tt,ss)/dot_product(tt,tt))
!11

      x_new=h_new+omega_new*ss
!12 if x accurate enough then quit


        if (all(abs(x_new-x_old)<=x_new*1.0E-10_dp)) then
          soln=x_new
          return
        else

          r_new=ss-omega_new*tt
        end if

    x_old=x_new
    rho_old=rho_new
    r_old=r_new
    omega_old=omega_new
    p_old=p_new
    v_old=v_new
    h_old=h_new
    !print*,x_new
end do
print*,'did not converge'



  end subroutine
end module
