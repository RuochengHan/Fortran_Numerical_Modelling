#define N 51

module deriv
contains
	function derivative_2(y, dx)
		implicit none
		real, intent(in)      :: y(N), dx
		real                  :: derivative_2(N)
  	integer               :: i
  	do i = 1, (N-1)
  		derivative_2(i) = (y(i+1) + y(i-1) - 2 * y(i)) / (dx * dx)
  	end do
	end function derivative_2
end module deriv

subroutine renew(k, dt, T, temp, dx)
  use deriv

  implicit none
  real                    :: ddy(N)
  real, intent(inout)     :: temp(N)
  real, intent(in)        :: k, dx, dt, T(N)
  integer                 :: i
	ddy = derivative_2(T, dx)
  do i = 2, N-1
    temp(i) = T(i) + dt * k * ddy(i) * dx * dx
  end do
end subroutine renew

program diff_1d
  implicit none
  integer                 :: i, choice
  real, parameter         :: k = 1.0
  real                    :: L, time, dx, total_time, dt
  real                    :: T(N), temp(N)
  ! ask and initialize the parameters
  print*, 'input length and total time:'
  read*, L, total_time
  time = 0.0
  dx = L / (N - 1)
  dt = 0.4 * (dx ** 2)

  ! initialize the T matrix
  call srand(94568)
  print*, 'choose initialization by random noise (1) or delta function (2):'
  read*, choice
	do i = 1, N
		T(i) = 0.0
	end do
	select case (choice)
	case (1)
		do i = 2, N-1
			T(i) = rand()
		end do
		print*, 'initialize with random noise'
	case (2)
		T(int((N+1)/2)) = 1.0
		print*, 'initialize with delta function'
	case default
		stop 'wrong input! type 1 or 2'
	end select

	! write initialized data to file
	open(1, file = 'diff_init.dat')
	do i = 1, N
		write(1, *) T(i)
	end do
	close(1)

  ! T matrix diffusion
  do while (time < total_time)
    call renew(k, dt, T, temp, dx)
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
end program diff_1d
