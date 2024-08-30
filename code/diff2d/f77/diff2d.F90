Program diff2d
  implicit none

#define REAL real*4

  integer :: m=4,n=4,b=1
  REAL,dimension(:,:),allocatable :: X,Y
  REAL :: norm
  
  allocate( X(-b:m+b-1,-b:n+b-1), Y(-b:m+b-1,-b:n+b-1) )

  block
    integer :: it,itcount
    do it=0,itcount-1
       call central_difference_from(x,y,m,n,b)
       norm = l2norm(x,m,n,b)
       call scale_interior(x,y,norm,m,n,b)
    end do
  end block

contains

  subroutine central_difference_from(y,x, m,n,b)
    implicit none
    ! paramteres
    integer,intent(in) :: m,n,b
    REAL,dimension(-b:m+b-1,-b:n+b-1),intent(in) :: x
    REAL,dimension(-b:m+b-1,-b:n+b-1),intent(out) :: y
    ! local variables
    integer :: row,col

    do row=0,m-1
       do col=0,n-1
          y(row,col) = 4*x(row,col) &
               - x(row+1,col) - x(row-1,col) &
               - x(row,col+1) - x(row,col-1)
       end do
    end do

  end subroutine central_difference_from

  function l2norm(x, m,n,b) result(norm)
    implicit none
    ! paramteres
    integer,intent(in) :: m,n,b
    REAL,dimension(-b:m+b-1,-b:n+b-1),intent(in) :: x
    REAL :: norm
    ! local variables
    integer :: row,col

    do row=0,m-1
       do col=0,n-1
          norm = norm + x(row,col)*x(row,col)
       end do
    end do

  end function l2norm

  subroutine scale_interior(y,x,norm, m,n,b)
    implicit none
    ! paramteres
    integer,intent(in) :: m,n,b
    REAL, intent(in) :: norm
    REAL,dimension(-b:m+b-1,-b:n+b-1),intent(in) :: x
    REAL,dimension(-b:m+b-1,-b:n+b-1),intent(out) :: y
    ! local variables
    integer :: row,col

    do row=0,m-1
       do col=0,n-1
          y(row,col) = x(row,col)/norm
       end do
    end do

  end subroutine scale_interior

End Program diff2d
