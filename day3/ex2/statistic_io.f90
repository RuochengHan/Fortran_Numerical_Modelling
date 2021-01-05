module stat
contains
	subroutine mean_stddev(arr, n, mean, stddev)
		implicit none
		real :: arr(n)
		real :: mean, stddev
		integer :: n, i
		real :: S, S_2
		i = 1
		S = 0.0
		S_2 = 0.0
		do while (i <= n)
			S = S + arr(i)
			S_2 = S_2 + arr(i) * arr(i)
			i = i + 1
		end do
		mean = S / n
		stddev = (S_2 / n - mean * mean)**(0.5)
  end subroutine mean_stddev
end module stat

program statistic_io
	use stat

	implicit none
	integer :: n, i
  real :: num, mean, stddev
  real, allocatable :: arr(:)
  i = 1
	n = 0
  mean = 0.0
  stddev = 0.0

	open(1, file = 'data.dat', status = 'old')
	do
  	read(1, *, iostat=i)
  	if (i < 0) exit
		if (i /= 0) stop 'error reading data'
  	n = n + 1
	end do
  allocate(arr(n))
	rewind(1)
	do i = 1, n
		read(1, *) arr(i)
	end do
	close(1)

  call mean_stddev(arr, n, mean, stddev)
	print*, "mean value is: ", mean, "standard deviation is ", stddev
  deallocate(arr)
end program statistic_io
