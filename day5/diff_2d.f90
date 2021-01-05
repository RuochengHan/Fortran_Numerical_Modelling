#define N 51
#define K 1
#define PI 3.1415926535

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
		del_squared(1, :) = 0.0
		del_squared(N, :) = 0.0

	end function del_squared

	function deriv_2d(T, vx, vy, h)
		implicit none
		real, intent(in)      :: T(N, N), vx(N, N), vy(N, N), h
		real                  :: deriv_2d(N, N)
		integer               :: i, j
  	do i = 2, N-1
			do j = 2, N-1
				deriv_2d(i, j) = 0.0
				if (vx(i, j) > 0) then
					deriv_2d(i, j) = deriv_2d(i, j) + vx(i, j) * (T(i, j) - T(i, j-1)) / h
				else
					deriv_2d(i, j) = deriv_2d(i, j) + vx(i, j) * (T(i, j+1) - T(i, j)) / h
				end if
				if (vy(i, j) > 0) then
					deriv_2d(i, j) = deriv_2d(i, j) + vy(i, j) * (T(i, j) - T(i-1, j)) / h
				else
					deriv_2d(i, j) = deriv_2d(i, j) + vy(i, j) * (T(i+1, j) - T(i, j)) / h
				end if
  		end do
		end do
	end function deriv_2d

	subroutine stream(S, vx, vy, h, vmax)
		implicit none
		real, intent(in)      :: S(N, N), h
		real, intent(inout)   :: vx(N, N), vy(N, N), vmax
		integer               :: i, j
		vmax = 0.0
		do i = 2, N-1
			do j = 2, N-1
    		vx(i, j) = (S(i+1, j) - S(i-1, j)) / (2 * h)
				vy(i, j) = - (S(i, j+1) - S(i, j-1)) / (2 * h)
				if (vx(i, j) > vmax) then
					vmax = vx(i, j)
				else if (vy(i, j) > vmax) then
					vmax = vy(i, j)
				end if
  		end do
		end do
	end subroutine stream
end module deriv

subroutine renew(T, v_dT, vx, vy, temp, h, dt)
  use deriv

  implicit none
  real                    :: ddT(N, N), v_dT(N, N)
  real, intent(inout)     :: temp(N, N)
  real, intent(in)        :: dt, h, T(N, N), vx(N, N), vy(N, N)
  integer                 :: i, j
	ddT = del_squared(T, h)
	v_dT = deriv_2d(T, vx, vy, h)
  do i = 1, N
		do j = 2, N-1
    temp(i, j) = T(i, j) + dt * (K * ddT(i, j) - v_dT(i, j))
		end do
  end do
	temp(:, 1) = temp(:, 2)
	temp(:, N) = temp(:, N-1)
end subroutine renew

program diff_2d
 use deriv

  implicit none
  integer                 :: i, j
  real                    :: L, a_diff, a_adv, B, h, time, total_time, dt, T(N, N), &
	                           temp(N, N), S(N, N), v_dT(N, N), vx(N, N), vy(N, N), vmax
	character(len=50)       :: choice, outputfile_1, outputfile_2

	! initialize the parameters from input file
	namelist /inputs/ L, a_diff, a_adv, B, total_time, choice, outputfile_1, outputfile_2
	open(3, file='params.dat')
	read(3, inputs)
	close(3)
	write(*, inputs)
  time = 0.0
  h = L / (N - 1)

	! initialize the S matrix
	do i = 1, N
		do j = 1, N
			S(i, j) = B * sin(PI * (i - 1) * h / L) * sin(PI * (j - 1) * h / L)
		end do
	end do

	! initialize the vx, vy matrix
	call stream(S, vx, vy, h, vmax)

	! initialize dt
	dt = min(a_diff * h * h / K, a_adv * h / vmax)

  ! initialize the T matrix
  call srand(94568)
	T(:, :) = 0.0
  select case (choice)
  case ('random')
		do i = 2, N-1
			do j = 1, N
      	T(i, j) = rand()
			end do
		end do
    print*, 'initialize with random noise'
  case ('delta')
		T(int((N+1)/2), int((N+1)/2)) = 1.0
    print*, 'initialize with delta function'
  case default
    stop 'wrong input!'
  end select
	T(1, :) = 1.0
	T(N, :) = 0.0

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
    call renew(T, v_dT, vx, vy, temp, h, dt)
		T = temp
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
