#define N 51
#define K 1

module deriv
contains
	function del_squared(T, h)
		implicit none
		real, intent(in)      :: T(N, N), h
		real                  :: del_squared(N, N)
		integer               :: i, j
  	do i = 2, (N-1)
			do j = 2, (N-1)
    		del_squared(i, j) = (T(i+1, j) + T(i-1, j) + T(i, j+1) + T(i, j-1) - 4 * T(i, j)) / (h * h)
  		end do
		end do
	end function del_squared
end module deriv

subroutine renew(T, temp, h, dt)
  use deriv

  implicit none
  real                    :: ddT(N, N)
  real, intent(inout)     :: temp(N, N)
  real, intent(in)        :: dt, h, T(N, N)
  integer                 :: i, j
	ddT = del_squared(T, h)
  do i = 2, N-1
		do j = 2, N-1
    temp(i, j) = T(i, j) + dt * K * ddT(i, j) * h * h
		end do
  end do
end subroutine renew

program diff_2d
  implicit none
  integer                 :: i, j
  real                    :: L, a, h, time, total_time, dt, T(N,N), temp(N,N)
	character(len=50)       :: choice, outputfile_1, outputfile_2

	! initialize the parameters from input file
	namelist /inputs/ L, a, total_time, choice, outputfile_1, outputfile_2
	open(3, file='params.dat')
	read(3, inputs)
	close(3)
	write(*, inputs)
  time = 0.0
  h = L / (N - 1)
  dt = a * (h ** 2) / K

  ! initialize the T matrix
  call srand(94568)
	do i = 1, N
		do j = 1, N
			T(i, j) = 0.0
		end do
	end do
  select case (choice)
  case ('random')
		do i = 2, N-1
			do j = 2, N-1
      	T(i, j) = rand()
			end do
		end do
    print*, 'initialize with random noise'
  case ('delta')
		T(int((N+1)/2), int((N+1)/2)) = 1.0
    print*, 'initialize with delta function'
  case default
    stop 'wrong input! type 1 or 2'
  end select

	! write initialized data to file
	open(1, file = outputfile_1)
	do i = 1, N
		do j = 1, N
			write(1, '(f10.6)', advance='no') T(i, j)
		end do
		write(1, *) ''
	end do
	close(1)

  ! T matrix diffusion
  do while (time < total_time)
    call renew(T, temp, h, dt)
    do i = 2, N-1
			do j = 2, N-1
      	T(i, j) = temp(i, j)
			end do
    end do
    time = time + dt
  end do

  ! write the result to file
  open(2, file = outputfile_2)
  do i = 1, N
		do j = 1, N
    	write(2, '(f10.6)', advance='no') T(i, j)
		end do
		write(2, *) ''
  end do
  close(2)
end program diff_2d
