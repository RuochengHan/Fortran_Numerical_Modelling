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
    real(kind=dp)                        :: h2, res_rms, sum
    integer                              :: i, j, nx, ny

    nx = size(u, 1)
    ny = size(u, 2)
    h2 = h * h
    sum = 0.0

    do i = 2, (nx-1)
      do j = 2, (ny-1)
         res_rms = (u(i+1, j) + u(i-1, j) + u(i, j+1) + u(i, j-1) - 4 * u(i, j)) / h2 - f(i, j)
         u(i, j) = u(i, j) + alpha * res_rms * h2 / 4
         sum = sum + res_rms * res_rms
      end do
    end do
    res_rms = sqrt(sum / ((nx - 2) * (ny - 2)))
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

module deriv
contains

  function deriv_2d_y(T, h) result (Tdy)
    implicit none
    integer, parameter                   :: dp = selected_real_kind(15, 307)
    real(kind=dp), intent(in)            :: h
    real(kind=dp), allocatable           :: T(:, :), Tdy(:, :)
    integer                              :: i, j, nx, ny
    nx = size(T, 1)
    ny = size(T, 2)
    allocate(Tdy(nx, ny))

    Tdy(:, :) = 0.0
    do i = 1, nx
      do j = 2, ny-1
         Tdy(i, j) = (T(i, j+1) - T(i, j-1)) / (2 * h)
      end do
    end do
  end function deriv_2d_y

	function del_squared(T, h)
		implicit none
    integer, parameter                   :: dp = selected_real_kind(15, 307)
		real(kind=dp), intent(in)            :: h
		real(kind=dp), allocatable           :: T(:, :), del_squared(:, :)
		integer                              :: i, j, nx, ny
		nx = size(T, 1)
		ny = size(T, 2)
		allocate(del_squared(nx, ny))

  	do i = 2, (nx-1)
			do j = 2, (ny-1)
    		del_squared(i, j) = (T(i+1, j) + T(i-1, j) + T(i, j+1) + T(i, j-1) - 4 * T(i, j)) / (h * h)
  		end do
		end do
		del_squared(1, :) = 0.0
		del_squared(nx, :) = 0.0
	end function del_squared

	function deriv_2d(T, vx, vy, h)
		implicit none
    integer, parameter                   :: dp = selected_real_kind(15, 307)
		real(kind=dp), intent(in)            :: h
		real(kind=dp), allocatable           :: T(:, :), vx(:, :), vy(:, :), deriv_2d(:, :)
		integer                              :: i, j, nx, ny
    nx = size(T, 1)
		ny = size(T, 2)
		allocate(deriv_2d(nx, ny))

  	do i = 2, nx-1
			do j = 2, ny-1
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
    integer, parameter                   :: dp = selected_real_kind(15, 307)
		real(kind=dp), intent(in)            :: h
		real(kind=dp), intent(inout)         :: vmax
		real(kind=dp), allocatable           :: S(:, :), vx(:, :), vy(:, :)
		integer                              :: i, j, nx, ny
    nx = size(S, 1)
		ny = size(S, 2)

		vmax = 0.0
		do i = 2, nx-1
			do j = 2, ny-1
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

program convect_2d
  use deriv
  use Poisson_Solver

  implicit none
  integer, parameter           :: dp = selected_real_kind(15, 307)
  integer                      :: i, j, k, Nx, Ny
  real(kind=dp)                :: Ra, a_diff, a_adv, kp, Pr, alpha, h, time, total_time, &
                                  w_rms, res, resall, thr, dt, vmax
  real(kind=dp), allocatable   :: T(:, :), tmp(:, :), Tdy(:, :), rhs(:, :), w(:, :), &
                                  res_mat(:, :), S(:, :), v_dT(:, :), vx(:, :), vy(:, :)
  character(len=50)            :: choice, outputfile_1, outputfile_2, outputfile_S, &
                                  outputfile_w, outputfile_T

	! initialize the parameters from input file
  namelist /inputs/ Nx, Ny, Ra, Pr, thr, kp, alpha, a_diff, a_adv, total_time, &
                    choice, outputfile_1, outputfile_2, outputfile_T, &
                    outputfile_w, outputfile_S
  open(3, file='params.dat')
	read(3, inputs)
	close(3)
	write(*, inputs)
  time = 0.0
  h = 1.0 / (Nx - 1)
  allocate(T(Nx, Ny), tmp(Nx, Ny), Tdy(Nx, Ny), rhs(Nx, Ny), w(Nx, Ny), S(Nx, Ny), &
           v_dT(Nx, Ny), vx(Nx, Ny), vy(Nx, Ny))

  ! initialize Temperature field (T matrix) with random or spike
  print*, "//////////////////////////////////////////////////"
  print*, "Initialize Temperature matrix"
  T(:, :) = 0.0
  select case (choice)
  case ('random')
    call random_number(T)
    print*, 'Initialize with random noise'
  case ('spike')
    T((Nx+1)/2, (Ny+1)/2) = 1.0 * Nx * Ny
    print*, 'Initialize with delta function'
  case default
    stop 'wrong input!'
  end select
  T(1, :) = 1.0
  T(Nx, :) = 0.0

  print*, "Initialize w matrix"
  w(:, :) = 0.0
  call random_number(w)
  print*, 'Initialize with random noise'
  print*, "//////////////////////////////////////////////////"
  write(*, *)

  ! write initialized data to file
  print*, "//////////////////////////////////////////////////"
  print*, "Write Initial Temperature matrix to file"
  open(1, file = outputfile_1)
  do i = 1, Nx
    do j = 1, Ny
      write(1, '(f15.6)', advance='no') T(i, j)
    end do
    write(1, *) ''
  end do
  close(1)

  print*, "Write Initial w matrix to file"
  open(2, file = outputfile_2)
  do i = 1, Nx
    do j = 1, Ny
      write(2, '(f15.6)', advance='no') w(i, j)
    end do
    write(2, *) ''
  end do
  close(2)
  print*, "//////////////////////////////////////////////////"
  write(*, *)

  ! Time Step
  print*, "//////////////////////////////////////////////////"
  !print*, "Start Time Step ..."
  do while (time < total_time)

    ! calculate rms of w
    w_rms = 0.0
    do i = 1, Nx
      do j = 1, Ny
        w_rms = w_rms + w(i, j) * w(i, j)
      end do
    end do
    w_rms = sqrt(w_rms / (Nx * Ny))

    ! multigrid possion: solve w for S
    res = 100.0
    k = 1
    do while(res > thr)
      res = abs(Vcycle_2DPoisson(S, w, h, alpha)) / w_rms
      !print*, k, " 's iteration residue is: ", res
      if (res > 1e5) then
        stop "Diverge ! Please try smaller alpha"
      end if
      k = k + 1
    end do

    ! calculate vx, vy from S
    call stream(S, vx, vy, h, vmax)

    ! calculate Tdy of T matrix and then rhs
    Tdy = deriv_2d_y(T, h)
    rhs = Ra * Tdy

    ! calculate dt under stable condition
    dt = min(a_diff * h * h / max(kp, kp * Pr), a_adv * h / vmax)

    ! diffusion of w field
    call renew(w, v_dT, vx, vy, tmp, h, kp, Pr, rhs, dt)

    ! diffusion of T field
    rhs = rhs * 0.0
    call renew(T, v_dT, vx, vy, tmp, h, kp, Pr / Pr, rhs, dt)

    ! show simulation time step
    time = time + dt
    print*, "simulation time is: ", time

  end do
  print*, "//////////////////////////////////////////////////"
  write(*, *)

  ! write last Stream matrix to file
  print*, "//////////////////////////////////////////////////"
  print*, "Write Final Stream matrix to file"
  open(3, file = outputfile_S)
  do i = 1, Nx
    do j = 1, Ny
      write(3, '(f15.6)', advance='no') S(i, j)
    end do
    write(3, *) ''
  end do
  close(3)
  print*, "//////////////////////////////////////////////////"
  write(*, *)

  ! write last Stream matrix to file
  print*, "//////////////////////////////////////////////////"
  print*, "Write Final w matrix to file"
  open(4, file = outputfile_w)
  do i = 1, Nx
    do j = 1, Ny
      write(4, '(f15.6)', advance='no') w(i, j)
    end do
    write(4, *) ''
  end do
  close(4)
  print*, "//////////////////////////////////////////////////"
  write(*, *)

  ! write last Temperature matrix to file
  print*, "//////////////////////////////////////////////////"
  print*, "Write Final Temperature matrix to file"
  open(5, file = outputfile_T)
  do i = 1, Nx
    do j = 1, Ny
      write(5, '(f15.6)', advance='no') T(i, j)
    end do
    write(5, *) ''
  end do
  close(5)
  print*, "//////////////////////////////////////////////////"

  deallocate(T, Tdy, rhs, w, S, v_dT, tmp, vx, vy)

contains
	subroutine renew(T, v_dT, vx, vy, tmp, h, kp, Pr, rhs, dt)
	  use deriv

	  implicit none
    integer, parameter                   :: dp = selected_real_kind(15, 307)
	  real(kind=dp), allocatable           :: ddT(:, :), v_dT(:, :), tmp(:, :), &
                                            T(:, :), vx(:, :), vy(:, :), rhs(:, :)
	  real(kind=dp), intent(in)            :: dt, h, kp, Pr
	  integer                              :: i, j, nx, ny
		nx = size(T, 1)
		ny = size(T, 2)
		allocate(ddT(nx, ny))

		ddT = del_squared(T, h)
		v_dT = deriv_2d(T, vx, vy, h)
	  do i = 1, nx
			do j = 2, ny-1
	    tmp(i, j) = T(i, j) + dt * (kp * Pr * ddT(i, j) - v_dT(i, j) - Pr * rhs(i, j))
			end do
	  end do
		tmp(:, 1) = tmp(:, 2)
		tmp(:, ny) = tmp(:, ny-1)
    T = tmp
		deallocate(ddT)
	end subroutine renew
end program convect_2d
