program neutronprog
  use constants
  use matrixmain
  use neutronsolve
  use bicg
  implicit none
  integer :: matrixsize
  type(matrix) :: matrixtest
  integer,dimension(:),allocatable :: row
  integer,dimension(:),allocatable :: col
  real(kind=dp),dimension(:),allocatable :: values

!for individual allocation ==================
  ! integer :: row
  ! integer :: col
  ! real(kind=dp) :: values
  !==========================================

  integer :: rsel,csel !to find a certain value
  !remember to delcare all variables
  integer :: ii
  integer,dimension(:),allocatable :: matrixrow
  real(kind=dp),dimension(:),allocatable :: multiplier
  real(kind=dp),dimension(:),allocatable :: result
  real(kind=dp),dimension(:),allocatable :: soln
    ! matrixsize = 5
    ! values=(/1,2,3,4,5,11,12,100,101/)
    ! row = (/1,2,3,4,5,1,2,4,5/)
    ! col = (/1,2,3,4,5,3,4,1,2/)

    ! matrixsize = 3
    ! values=(/1,101,2,100,3/)
    ! row = (/1,1,2,3,3/)
    ! col = (/1,3,2,1,3/)

    matrixsize = 7
    allocate(values(matrixsize))
    allocate(row(matrixsize))
    allocate(col(matrixsize))
    values=(/1,2,3,4,5,6,7,-10,-11,-12,100,1000,2000,-100,-101,-102,4,1,1/)
    row = (/1,2,3,4,5,6,7,1,2,3,1,2,3,2,3,4,7,6,7/)
    col = (/1,2,3,4,5,6,7,3,4,5,5,6,7,1,2,3,6,2,3/)


    !print*,values,row,col

  call matrixtest%initialise(matrixsize)
do ii=1,size(values)
  call matrixtest%addvalue(row(ii),col(ii),values(ii))

!for adding in loop========================
! do ii=1,100
!   row=ii
!   col=ii
!   values=ii
!
!   call matrixtest%addvalue(row,col,values)
!==========================================
end do

!print*, matrixtest%getvalue(2,1)


call matrixtest%build()
allocate(multiplier(matrixsize))
allocate(result(matrixsize))
!multiplier=(/10.0_dp,11.0_dp,12.0_dp,13.0_dp,15.0_dp/)
multiplier=(/10.0_dp,11.0_dp,12.0_dp,13.0_dp,14.0_dp,15.0_dp,16.0_dp/)

!multiplier=(/1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp/)
!print*,multiplier
call matrixtest%multiply(multiplier,result)
!print*, result

!print*,'idx',matrixtest%idx
!print*,'values',matrixtest%values
!print*,'ndiag',matrixtest%numdiag
allocate(soln(matrixsize))
!soln=(/0.99647023027119752_dp,2.0001300821105583_dp,2.9985425543148110_dp,3.9908538604945476_dp,5.0004140951727001_dp/)
call bicgsolve(matrixtest,result,soln)


print*,'soln', soln
end program
