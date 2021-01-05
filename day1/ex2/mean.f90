program mean
	implicit none
	real num1, num2, num3
	print*, "input three numbers"
	read*, num1
	read*, num2
	read*, num3

	!calculate arithmetic mean, harmonic mean and geometric mean
	print*, "arithmetic mean is ", (num1 + num2 + num3)/3.0
	print*, "harmonic mean is ", 3.0/(1.0/num1 + 1.0/num2 + 1.0/num3)
	print*, "geometric mean is ", (num1 * num2 * num3)**(1.0/3.0)
end program mean

