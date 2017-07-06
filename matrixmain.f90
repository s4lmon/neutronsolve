
! line 92 not allocating correclty
!
!



module matrixmain
  use constants
  implicit none


  type :: matrix
    real(kind=dp),dimension(:,:),allocatable :: values
    integer :: numdiag
    integer,dimension(:),allocatable :: idx
    integer :: matrixsize
  contains
    procedure :: initialise=>initialise_matrix
    procedure :: addvalue=>addvalue_matrix
    procedure :: getvalue=>getvalue_matrix
    procedure :: build=>build_matrix
    procedure :: multiply=>multiply_matrix
    !procedure :: valueset=>valuematrix
    !procedure :: initialise=>initialisesquare
  end type

contains
  subroutine initialise_matrix(this,matrixsize)
    use constants
    class(matrix),intent(inout) :: this
    integer :: matrixsize !number of rows = number of columns
  !starting with it at zero
    this%numdiag=0
    !setting to 0 but will always be same as numdiag
    allocate(this%idx(this%numdiag))
    !creating blank future compressed diagonal matrix
    allocate(this%values(matrixsize,this%numdiag))
    this%matrixsize = matrixsize
  end subroutine

  subroutine addvalue_matrix(this,row,col,values)
    use constants
    class(matrix),intent(inout) :: this
    integer :: col
    integer :: row
    integer :: new_diag !this is finding the index (e.g. all row=column values are on leading diagonal)
    integer :: ii
    integer,dimension(:),allocatable :: cidx !
    real(kind=dp),dimension(:,:),allocatable :: cvalues !test
    integer :: loc
    real(kind=dp),dimension(:,:),allocatable :: csm !not sure if needed
    real(kind=dp) :: values
    integer :: col_num !used for adding to values matrix, natural number unlike index
    !xx!xxprint *, col,row,values
    new_diag = col-row

    if(any(this%idx==new_diag)) then !first if to check if any diags already exist
      !allocate(cidx(size(this%idx))) !not essential
      !cidx=this%idx not needed, in this case entry is already there
      !allocate(cvalues(this%matrixsize,size(this%idx)))

      col_num=minloc(this%idx,dim=1,mask=this%idx .eq. new_diag)
      !xx *, 'found','colnum',col_num
      this%values(row,col_num)=values

      !------------------------------
      !------------------------------
      else

      !xx!xxprint *, 'not found'
      allocate(cidx(size(this%idx)+1))
      allocate(cvalues(this%matrixsize,size(this%idx)+1))
      cvalues=0.0 !always allocate reals with .0
      !------------------------------
      !when no value is found
          if (size(this%idx)==0) then !if matrix is already empty, special case
              cidx(1)=new_diag !omitting (1) has no effect because matrix already size 1
              this%numdiag=1
              col_num=minloc(cidx,dim=1,mask=cidx .eq. new_diag)
              cvalues(row,col_num)=values
              this%numdiag=this%numdiag+1 !fix
              !this%values(row,cidx)=values

                          !this%values(1)=values(ii)
          else
            if  (minloc(this%idx,dim=1,mask=this%idx > new_diag)==0) then
              col_num=size(this%idx)+1
              this%numdiag=this%numdiag+1 !fix
            else
            col_num=minloc(this%idx,dim=1,mask=this%idx > new_diag) ! should be >
          end if
            !xxprint*, 'test',this%idx,new_diag, col_num
            cidx(:col_num-1)=this%idx(:col_num-1)
            cidx(col_num)=new_diag
            cidx(col_num+1:)=this%idx(col_num:)
            cvalues(:,:col_num-1)=this%values(:,:col_num-1)
            cvalues(row,col_num)=values
            cvalues(:,col_num+1:)=this%values(:,col_num:)






            !old logic
              ! do ii=1,size(cidx)
              !   if (new_diag < this%idx(ii)) then
              !     cidx(ii)=new_diag
              !
              !     cvalues(row,col_num)=values
              !       !xxprint *, cidx(ii)
              !     cidx(ii+1:)=this%idx(ii:) !taking values after place where new diag fits in
              !     !this%values(row,cidx)=values(ii)
              !     exit
              !   else
              !     cidx(ii)=this%idx(ii)
              !     cvalues(row,col_num)=values !check this
              !
              !     loc = ii
              !   end if
              ! end do
          end if
              deallocate(this%idx)
              allocate(this%idx(size(cidx)))
              this%idx=cidx
              deallocate(this%values)
              allocate(this%values(this%matrixsize,size(cvalues)))
              this%values=cvalues
              !------------------------------

              !------------------------------

    end if
    !xxprint *, 'idx',this%idx
    !xxprint *, 'values',this%values
      !now want size of values relayed back to main function to see if matrix input lies in this

  end subroutine
  ! subroutine valueset(this,)
  !   real(kind=dp)
  !   use constants
  !   do i=1,size(values)
  !     call matrix%addvalue
  !   end do
  !
  !   class(matrix)
  ! end subroutine


  ! function areasquare(this) result(evaluation)
  !
  !   class(square),intent(in) :: this
  !   real(kind=dp) :: evaluation
  !   evaluation=this%dimensions(1)*this%dimensions(2)
  ! end function
  !
  ! function perimetersquare(this) result(evaluation)
  !
  !   class(square),intent(in) :: this
  !   real(kind=dp) :: evaluation
  !
  ! end function


  function getvalue_matrix(this,rsel,csel) result(found)
    !do i need intent?
    class(matrix), intent(in) :: this
    real(kind=dp) :: found
    integer :: rsel,csel


    if(any(this%idx==(csel-rsel))) then

      found=this%values(rsel,(minloc(this%idx,dim=1,mask=this%idx == (csel-rsel))))

    else
      found=0
    end if
  end function

  subroutine build_matrix(this)

    class(matrix), intent(in) :: this

    integer :: max
    integer :: count
    integer :: count2
    real(kind=dp),dimension(:),allocatable :: matrixrow


    allocate(matrixrow(this%matrixsize))
    max=101!what's max matrix size to !xxprint?

    if (this%matrixsize>max) then
      print*, 'matrix is too large'
    else
      do count=1,this%matrixsize
        do count2=1,this%matrixsize
        matrixrow(count2)=this%getvalue(count,count2)

        end do
        print*,'row',count,'  ',matrixrow
      end do
    end if
  end subroutine

  subroutine multiply_matrix(this,multiplier,result)
    class(matrix), intent(in) :: this
    real(kind=dp),dimension(this%matrixsize) :: multiplier
    real(kind=dp),dimension(this%matrixsize) :: result
    integer :: cc
    integer :: rr

    result = 0


    do cc=1,size(this%idx)
      !write(*,'(x,a,i0)') "cc = ", cc
      do rr=max(1,1-this%idx(cc)),min(this%matrixsize,this%matrixsize-this%idx(cc))
      result(rr) = result(rr) + this%values(rr,cc)*multiplier(this%idx(cc)+rr)
      !write(*,'(x,a,3(i0,a))') "Multiplying row ", rr, " column ", cc, " of values by row ", this%idx(cc)+rr, " of multiplier."
    end do

    end do

  end subroutine


end module
