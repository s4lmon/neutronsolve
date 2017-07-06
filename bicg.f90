module bicg
  use constants
  use matrixmain
  implicit none

contains
  subroutine bicgsolve(matrixin,rhs,soln,iter)
    class(matrix),intent(inout) :: matrixin !no longer class-bound subroutine
    real(kind=dp),dimension(matrixin%matrixsize) :: r_0,rhs,soln,r_old,r_new,r_hat,v_new,v_old,x_old,x_new,p_old,p_new, h_old, h_new
    real(kind=dp),dimension(matrixin%matrixsize) :: ss,tt
    real(kind=dp) :: err !can't give intent to variables not in argumentlist
    real(kind=dp) :: alpha,rho_old,rho_new, beta, omega_old,omega_new
    integer :: count,iter,maxiter !test
    !logical :: converged
    v_old=0.0_dp
    p_old=0.0_dp
    soln=1.0_dp !test only, soln set arbitrarily from meror
    call matrixin%multiply(soln,r_0)
    r_0=rhs-r_0
    r_old=r_0
    r_hat=r_0
    rho_old=1.0_dp
    omega_old=1.0_dp
    alpha=1.0_dp
    x_old=soln

    maxiter=10000
do count=1,maxiter
  !print*, count, "===================="
!1

    rho_new=dot_product(r_hat,r_old)

!2

    beta=(rho_new/rho_old)*(alpha/omega_old)

!3

    p_new=r_old+beta*(p_old-omega_old*v_old)
!4
    call matrixin%multiply(p_new,v_new)


!5  !issue lies with p_new
    alpha=rho_new/(dot_product(r_hat,v_new))


!6

    h_new = x_old + (alpha*p_new)

!7
    if (all(abs(h_new-h_old)<=abs(h_new*1.0E-5_dp))) then !h is accurate enough
    ! converged = .true.
    ! do test=1, size(h_new)
    !   if (abs(h_new(test)-h_old(test)) .gt. h_new(test)*1e-10_dp) converged = .false.
    !   write(*,'(x,2(a,es12.4))') "h_new - h_old = ", h_new(test)-h_old(test), "exit criterion = ", h_new(test)*1e-10_dp
    ! end do
    ! if (converged) then

      soln=h_new

      iter=count
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


        if (all(abs(x_new-x_old)<=abs(x_new*1.0E-5_dp))) then


        ! converged = .true.
        ! do test=1, size(x_new)
        !   if (abs(x_new(test)-x_old(test)) .gt. abs(x_new(test)*1e-10_dp)) converged = .false.
        !   write(*,'(x,2(a,es12.4))') "x_new - x_old = ", x_new(test)-x_old(test), "exit criterion = ", x_new(test)*1e-10_dp
        ! end do
        ! if (converged) then
          soln=x_new
          iter=count
          return
        else

          r_new=ss-(omega_new*tt)
        end if

    x_old=x_new
    rho_old=rho_new
    r_old=r_new
    omega_old=omega_new
    p_old=p_new
    v_old=v_new
    h_old=h_new

end do
print*,'did not converge'



  end subroutine
end module
