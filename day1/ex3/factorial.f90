program factorial
	implicit none
	integer(kind=8) i, num, F
	F = 1
	do while (num < 1)
		print *, "input a positive integer: "
		read*, num
		if (num < 1) print *, "number not positive!"
	end do
	do i = 1,num
		F = F * i
	end do
	print *, num, "'s factorial is ", F
end program factorial
	
