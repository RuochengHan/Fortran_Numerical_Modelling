module Poisson_Solver
contains

  recursive function Vcycle_2DPoisson(u_f,rhs,h,alpha) result (resV)
    implicit none
    integer, parameter         :: dp = selected_real_kind(15, 307)
    real(kind=dp) resV
    real(kind=dp),intent(inout):: u_f(:,:)  ! arguments
    real(kind=dp),intent(in)   :: rhs(:,:),h
    integer                    :: nx,ny,nxc,nyc, i,j  ! local variables
    real(kind=dp),allocatable  :: res_c(:,:),corr_c(:,:),res_f(:,:),corr_f(:,:)
    real(kind=dp)              :: alpha, res_rms

    nx=size(u_f,1); ny=size(u_f,2)  ! must be power of 2 plus 1
    nxc=1+(nx-1)/2; nyc=1+(ny-1)/2  ! coarse grid size

    if (min(nx,ny)>5) then  ! not the coarsest level

       allocate(res_f(nx,ny),corr_f(nx,ny), &
            corr_c(nxc,nyc),res_c(nxc,nyc))

       !---------- take 2 iterations on the fine grid--------------
       res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)
       res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)

       !---------- restrict the residue to the coarse grid --------
       call residue_2DPoisson(u_f,rhs,h,res_f)
       call restrict(res_f,res_c)

       !---------- solve for the coarse grid correction -----------
       corr_c = 0.
       res_rms = Vcycle_2DPoisson(corr_c,res_c,h*2,alpha) ! *RECURSIVE CALL*

       !---- prolongate (interpolate) the correction to the fine grid
       call prolongate(corr_c,corr_f)

       !---------- correct the fine-grid solution -----------------
       u_f = u_f - corr_f

       !---------- two more smoothing iterations on the fine grid---
       res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)
       res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)

       deallocate(res_f,corr_f,res_c,corr_c)

    else

       !----- coarsest level (ny=5): iterate to get 'exact' solution

       do i = 1,100
          res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)
       end do

    end if

    resV = res_rms   ! returns the rms. residue

  end function Vcycle_2DPoisson

  function iteration_2DPoisson(u, f, h, alpha) result (res_rms)
    implicit none
    integer, parameter                   :: dp = selected_real_kind(15, 307)
    real(kind=dp), intent(inout)         :: u(:, :)
    real(kind=dp), intent(in)            :: f(:, :), h, alpha
    real(kind=dp)                        :: h2, res_rms
    integer                              :: i, j, nx, ny

    nx = size(u, 1)
    ny = size(u, 2)
    h2 = h * h

    do i = 2, (nx-1)
      do j = 2, (ny-1)
         res_rms = (u(i+1, j) + u(i-1, j) + u(i, j+1) + u(i, j-1) - 4 * u(i, j)) / h2 - f(i, j)
         u(i, j) = u(i, j) + alpha * res_rms * h2 / 4
      end do
    end do
  end function iteration_2DPoisson

  subroutine residue_2DPoisson(u, f, h, res)
    implicit none
    integer, parameter                   :: dp = selected_real_kind(15, 307)
    real(kind=dp), intent(inout)         :: u(:, :), res(:, :)
    real(kind=dp), intent(in)            :: f(:, :), h
    real(kind=dp)                        :: h2
    integer                     :: i, j, nx, ny

    nx = size(u, 1)
    ny = size(u, 2)
    h2 = h * h

    do i = 2, (nx-1)
      do j = 2, (ny-1)
         res(i, j) = (u(i+1, j) + u(i-1, j) + u(i, j+1) + u(i, j-1) - 4 * u(i, j)) / h2 - f(i, j)
      end do
    end do
  end subroutine residue_2DPoisson

  subroutine restrict(fine, course)
    implicit none
    integer, parameter                   :: dp = selected_real_kind(15, 307)
    real(kind=dp), intent(inout)         :: course(:, :)
    real(kind=dp), intent(in)            :: fine(:, :)
    integer                              :: i, j, nxc, nyc

    nxc = size(course, 1)
    nyc = size(course, 2)

    do i = 1, nxc
      do j = 1, nyc
         course(i, j) = fine((i-1)*2+1, (j-1)*2+1)
      end do
    end do
  end subroutine restrict

  subroutine prolongate(course, fine)
    implicit none
    integer, parameter                   :: dp = selected_real_kind(15, 307)
    real(kind=dp), intent(in)            :: course(:, :)
    real(kind=dp), intent(inout)         :: fine(:, :)
    integer                              :: i, j, nxc, nyc, nyf

    nxc = size(course, 1)
    nyc = size(course, 2)
    nyf = size(fine, 2)

    do i = 1, nxc
      do j = 1, nyc
         fine((i-1)*2+1, (j-1)*2+1) = course(i, j)
      end do
    end do

    do i = 1, nxc
      do j = 2, nyc
         fine((i-1)*2+1, (j-1)*2) = (course(i, j-1) + course(i, j)) / 2
      end do
    end do

    do i = 2, nxc
      do j = 1, nyf
         fine((i-1)*2, j) = (fine((i-1)*2-1, j) + fine((i-1)*2+1, j)) / 2
      end do
    end do
  end subroutine prolongate
end module Poisson_Solver

program test_poisson
  use Poisson_Solver

  implicit none
  integer, parameter                     :: dp = selected_real_kind(15, 307)
  real(kind=dp), allocatable             :: u_f(:, :), rhs(:, :), resI_f(:, :)
  real(kind=dp)                          :: h, resV, resI, thr, rhs_rms, rms, alpha, L
  real                                   :: t1, t2
  integer                                :: i, j, k, N
  character(len=50)                      :: choice, option, outputfile_1, outputfile_2

  ! initialize the parameters from input file
  namelist /inputs/ N, L, thr, alpha, choice, option, outputfile_1, outputfile_2
  open(3, file='params.dat')
  read(3, inputs)
  close(3)
  write(*, inputs)
  allocate(u_f(N, N), rhs(N, N), resI_f(N, N))
  h = L / (N - 1)

  ! initialize u_f matrix with all zero
  u_f(:, :) = 0.0

  ! initialize field (rhs matrix) with random or spike
  rhs(:, :) = 0.0
  select case (choice)
  case ('random')
    call random_number(rhs)
    print*, 'initialize with random noise'
  case ('spike')
    rhs((N+1)/2, (N+1)/2) = 1.0 * N * N
    print*, 'initialize with delta function'
  case default
    stop 'wrong input!'
  end select
  rhs(1, :) = 0.0
  rhs(N, :) = 0.0
  rhs(:, 1) = 0.0
  rhs(:, N) = 0.0

  rhs_rms = 0.0
  do i = 1, N
    do j = 1, N
      rhs_rms = rhs_rms + rhs(i, j) * rhs(i, j)
    end do
  end do
  rhs_rms = sqrt(rhs_rms) / N

  ! write field (rhs matrix) to file
  open(1, file = outputfile_1)
  do i = 1, N
    do j = 1, N
      write(1, '(f10.6)', advance='no') rhs(i, j)
    end do
    write(1, *) ''
  end do
  close(1)

  ! iterate using Vcycle or Iteration until meet the converge threshold
  call cpu_time(t1)

  select case (option)
  case ('Vcycle')
    print*, "Start Vcycle !"
    rms = 100.0
    k = 1
    do while(rms > thr)
      resV = abs(Vcycle_2DPoisson(u_f, rhs, h, alpha))
      rms = resV / rhs_rms
      print*, k, " 's iteration residue is: ", rms
      if (rms > 1e8) then
        stop 'Diverge ! Please try smaller alpha'
      end if
      k = k + 1
    end do
    print*, "Meet the converge threshold !"
  case ('Itercycle')
    print*, "Start Itercycle !"
    rms = 100.0
    k = 1
    do while(rms > thr)
      resI = abs(iteration_2DPoisson(u_f, rhs, h, alpha))
      call residue_2DPoisson(u_f, rhs, h, resI_f)
      rms = 0.0
      do i = 1, N
        do j = 1, N
          rms = rms + abs(resI_f(i, j))
        end do
      end do
      rms = rms / rhs_rms
      print*, k, " 's iteration residue is: ", rms
      if (rms > 1e8) then
        stop 'Diverge ! Please try smaller alpha'
      end if
      k = k + 1
    end do
    print*, "Meet the converge threshold !"
  case default
    stop 'wrong input!'
  end select

  call cpu_time(t2)
  print*, "Iteration time is ", (t2 - t1) * 1000, " ms"

  ! write the result to file
  open(2, file = outputfile_2)
  do i = 1, N
    do j = 1, N
      write(2, '(f10.6)', advance='no') u_f(i, j)
    end do
    write(2, *) ''
  end do
  close(2)

  deallocate(u_f, rhs, resI_f)
end program test_poisson
