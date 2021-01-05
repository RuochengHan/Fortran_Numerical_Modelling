program statistic
	implicit none
	integer n, i
	real num, S, S_2
	num = 0.0
	S = 0.0
	S_2 = 0.0
	i = 0
	print*, "input number of data: "
	read*, n
	do while (i < n)
		print*, "input a positive real number: "
		read*, num
		if (num <= 0.0) then
			print*, "number not positive!"
		else
		S = S + num
		S_2 = S_2 + num * num
		i = i + 1
		end if
	end do
	print*, "mean value is: ", S / n, "standard deviation is ", (S_2 / n - (S / n) * (S / n))**(0.5)
end program statistic
