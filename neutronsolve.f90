!do i do it for general case of is line 26 ok starting no 2
!is this format (module, subroutien) correct?

module neutronsolve
  use constants
  use matrixmain
  use bicg

  implicit none

contains
  subroutine nte(length,disc,sigma_t,source,x_bar,psi_bar)
    real(kind=dp) :: length !length of domain (m)
    real(kind=dp) :: sigma_t !macroscopic total cross section, which includes all possible interactions
    real(kind=dp) :: source !source of flux
    real(kind=dp) :: delta !size of space
    real(kind=dp),dimension(:),allocatable :: x_bar,psi_bar
    integer :: disc   !number of discretisations of domain = matrixsize
    !integer :: matrixsize
    type(matrix) :: matrixsolve
    integer  :: row
    integer  :: col,ii

    real(kind=dp) :: value
    delta=length/disc
    allocate(x_bar(disc))
    allocate(psi_bar(disc))
    x_bar=0 !first value = 0

    call matrixsolve%initialise(disc)
    call matrixsolve%addvalue(1,1,((2/delta)-sigma_t))
      do ii = 2,disc !1 is special case
          x_bar(ii)=ii*delta
          row = ii
          col = ii
          value = (-1/delta)-sigma_t
        call matrixsolve%addvalue(row,col,value)
          row = ii
          col = ii-1
          value = 1/delta
        call matrixsolve%addvalue(row,col,value)
      end do
  end subroutine
end module
