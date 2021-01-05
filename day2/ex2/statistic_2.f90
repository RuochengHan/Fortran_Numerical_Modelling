program statistic
	implicit none
	integer :: n, i
  real :: num, mean, stddev
  real, allocatable :: arr(:)
  i = 1
  mean = 0.0
  stddev = 0.0
	print*, "input number of data: "
	read*, n
  allocate(arr(n))
	do while (i <= n)
		print*, "input a positive real number: "
		read*, num
		if (num <= 0.0) then
			print*, "number not positive!"
    else
      arr(i) = num
      i = i + 1
		end if
	end do
  call mean_stddev(arr, n, mean, stddev)
	print*, "mean value is: ", mean, "standard deviation is ", stddev
  deallocate(arr)
end program statistic

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
