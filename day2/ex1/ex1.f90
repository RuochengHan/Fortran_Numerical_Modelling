program ex1
  implicit none

  character(15) :: name

  integer, parameter :: a = 5

  real, dimension(-1:10) :: arr1

  real, allocatable :: arr2(:,:,:,:)

  integer, parameter :: b = nint(3.6)

  integer :: c = mod(a, b)

  integer :: i, j, sum = 0
  real, dimension(1:100) :: arr3
  do i = 12,124,2
    sum = sum + i
  end do

  call random_number(arr3)
  do j = 1,100
    if (arr3(j) > 0) then
      print*, 'input a positive num ', arr3(j), "at No. ", j
      exit
    end if
  end do
end program ex1
