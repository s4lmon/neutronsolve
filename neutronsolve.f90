!do i do it for general case of is line 26 ok starting no 2
!is this format (module, subroutien) correct?

module neutronsolve
  use constants
  use matrixmain
  use bicg

  implicit none

contains
  subroutine nte(length,disc,sigma_t,source,psi_array,iter,x_array)
    type(matrix) :: matrixsolve

    real(kind=dp) :: length !length of domain (m)
    real(kind=dp) :: sigma_t !macroscopic total cross section, which includes all possible interactions
    real(kind=dp) :: source !source of flux
    real(kind=dp) :: delta !size of space
    real(kind=dp),dimension(:),allocatable :: x_array,psi_array,rhs
    integer :: disc,iter   !number of discretisations of domain = matrixsize
    !integer :: matrixsize

    integer  :: row
    integer  :: col,ii

    real(kind=dp) :: value
    delta=length/(disc-1)
    allocate(x_array(disc))
    allocate(psi_array(disc))
    allocate(rhs(disc))
    x_array=0 !first value = 0
    psi_array=0.0_dp
    call matrixsolve%initialise(disc)
    call matrixsolve%addvalue(1,1,((2.0_dp/delta)-sigma_t))



      do ii = 2,disc !1 is special case
          x_array(ii)=(ii-1)*delta
          row = ii
          col = ii
          value = (-1/delta)-sigma_t
        call matrixsolve%addvalue(row,col,value)
          row = ii
          col = ii-1
          value = 1/delta
        call matrixsolve%addvalue(row,col,value)
      end do

    !defining solution matrix, rhs
    rhs=0
    rhs(1)=(2*source)/delta

    !print*,'x_array',x_array
    call matrixsolve%build()
    !print*,'rhs',rhs
    call bicgsolve(matrixsolve,rhs,psi_array,iter)
  end subroutine
end module
