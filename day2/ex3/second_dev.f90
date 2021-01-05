program deriv2
  implicit none
  real, allocatable     :: x(:), y(:), ddy(:)
  real                  :: dx
  integer               :: n, i
  print*, 'input num of grid points: '
  read*, n
  allocate(x(n), y(n), ddy(n))
  dx = 10.0 / (n)
  do i = 0 ,n
    x(i) = i * dx
    y(i) = x(i) * x(i)
	end do
  call derivative_2(y, n, dx, ddy)
  do i = 1, (n-1)
    print*, ddy(i), 2, ddy(i) - 2
	end do
  deallocate(x, y, ddy)
end program deriv2

subroutine derivative_2(y, n, dx, ddy)
  real, intent(in)      :: y(n), dx
  integer, intent(in)   :: n
  real, intent(out)     :: ddy(n)
  integer               :: i
  do i = 1, (n-1)
    ddy(i) = (y(i+1) + y(i-1) - 2 * y(i)) / (dx * dx)
  end do
end subroutine derivative_2
