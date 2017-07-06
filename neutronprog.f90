program neutronprog
  use neutronsolve

  implicit none
  real(kind=dp) :: lengthofdomain
  integer :: numdiscret,iter,ii
  real(kind=dp) :: sigma_t
  real(kind=dp) :: source
  real(kind=dp),dimension(:),allocatable :: x_array,psi_array,rhs
lengthofdomain=10
numdiscret=3000
sigma_t=1
source=100



call nte(lengthofdomain,numdiscret,sigma_t,source,psi_array,iter,x_array)

!sx_array,psi_array)
  if (numdiscret<5000) then
    print*,'psi_array:'
    do ii=1,size(psi_array)
      print*, psi_array(ii),x_array(ii)
    end do
    print*,iter,"iterations"
    !print*,x_array
    ! do ii = 1,size(x_array)
    !   print*, x_array(ii)
    ! end do
  else
    print*,'results too large to print'
  end if
end program
