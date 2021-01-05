module deriv
contains
	subroutine derivative_2(y, n, dx, ddy)
	implicit none
	real, intent(in)      :: y(n), dx
  	integer, intent(in)   :: n
	real, intent(out)     :: ddy(n)
  	integer               :: i
  	do i = 1, (n-1)
    	ddy(i) = (y(i+1) + y(i-1) - 2 * y(i)) / (dx * dx)
  	end do
	end subroutine derivative_2
end module deriv

subroutine renew(k, N, T, temp, dx)
  use deriv

  implicit none
  real                    :: ddy(N)
  real, intent(inout)     :: temp(N)
  real, intent(in)        :: k, dx, T(N)
  integer, intent(in)     :: N
  integer                 :: i
  call derivative_2(T, N, dx, ddy)
  do i = 2, N-1
    temp(i) = T(i) + k * ddy(i) * dx * dx
  end do
end subroutine renew

program diff_1d
  implicit none
  integer                 :: i, N, choice
  real, parameter         :: k = 1.0
  real                    :: L, time, dx, total_time, dt
  real, allocatable       :: T(:), temp(:)
  ! ask and initialize the parameters
  print*, 'input length, discrete point and total time:'
  read*, L, N, total_time
  time = 0.0
  dx = L / (N - 1)
  dt = 0.4 * (dx ** 2)
  allocate(T(N), temp(N))

  ! initialize the T matrix
  call srand(94568)
  print*, 'choose initialization by random noise (1) or delta function (2):'
  read*, choice
  T(1) = 0.0
  T(N) = 0.0
  temp(1) = 0.0
  temp(N) = 0.0
  open(1, file='diff_init.dat')

  select case (choice)
  case (1)
    write(1, *) T(1)
    do i = 2, N-1
      T(i) = rand()
      write(1, *) T(i)
    end do
    write(1, *) T(N)
    close(1)
    print*, 'initialize with random noise'
  case (2)
    write(1, *) T(1)
    do i = 2, N-1
      T(i) = exp(-(i * dx - L / 2) * (i * dx - L / 2) / 1e-3)
      write(1, *) T(i)
    end do
    write(1, *) T(N)
    close(1)
    print*, 'initialize with delta function'
  case default
    stop 'wrong input! type 1 or 2'
  end select

  ! T matrix diffusion
  do while (time < total_time)
    call renew(k, N, T, temp, dx)
    do i = 2, N-1
      T(i) = temp(i)
      !print*, T(i)
    end do
    time = time + dt
  end do

  ! output the result to file
  open(2, file='diff_result.dat')
  do i = 1, N
    write(2, *) T(i)
  end do
  close(2)
  deallocate(T, temp)
end program diff_1d
