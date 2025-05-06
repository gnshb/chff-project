module data
  implicit none
  integer, parameter :: nx = 60, ny = 60
  integer, parameter :: u_iter = 30, v_iter = 30, p_iter = 250
  double precision, parameter :: Lx = 1.0d0, Ly = 1.0d0
  double precision, parameter :: mu = 0.001d0, ulid = 1.0d0, rho = 1.0d0
  double precision, parameter :: alpha_u = 0.70d0, alpha_v = 0.70d0, alpha_p = 0.30d0
  double precision :: u_residual, v_residual, p_residual, residual 
  double precision :: dx, dy 
  double precision, dimension(:, :), allocatable :: uf, vf, p, xgrid, ygrid, uc, vc
  double precision, dimension(:, :), allocatable :: uf_old, vf_old, pprime
  double precision, dimension(:, :), allocatable :: ap_u, ae_u, aw_u, an_u, as_u, b_u, d_u
  double precision, dimension(:, :), allocatable :: ap_v, ae_v, aw_v, an_v, as_v, b_v, d_v
  double precision, dimension(:, :), allocatable :: ap_p, ae_p, aw_p, an_p, as_p, b_p
  
contains 
  
  subroutine allocate_memory()
    implicit none
    
    allocate (xgrid(0:nx+1, 0:ny+1), ygrid(0:nx+1, 0:ny+1))
    
    allocate (p(0:nx+1, 0:ny+1))
    allocate (uf(0:nx, 0:ny+1))
    allocate (vf(0:nx+1, 0:ny))
    
    allocate (uc(0:nx+1, 0:ny+1), vc(0:nx+1, 0:ny+1))

    allocate (pprime(0:nx+1, 0:ny+1))
    allocate (uf_old(0:nx, 0:ny+1))
    allocate (vf_old(0:nx+1, 0:ny))
    
    allocate (ap_p(1:nx, 1:ny), ae_p(1:nx, 1:ny), aw_p(1:nx, 1:ny), &
         &    an_p(1:nx, 1:ny), as_p(1:nx, 1:ny), b_p(1:nx, 1:ny))
    
    allocate (ap_u(0:nx, 1:ny), ae_u(0:nx, 1:ny), aw_u(0:nx, 1:ny), &
         &    an_u(0:nx, 1:ny), as_u(0:nx, 1:ny), b_u(0:nx, 1:ny), d_u(0:nx, 1:ny))
    
    allocate (ap_v(1:nx, 0:ny), ae_v(1:nx, 0:ny), aw_v(1:nx, 0:ny), &
         &    an_v(1:nx, 0:ny), as_v(1:nx, 0:ny), b_v(1:nx, 0:ny), d_v(1:nx, 0:ny))
    
  end subroutine allocate_memory
  
  subroutine deallocate_memory()
    implicit none
    
    deallocate (xgrid, ygrid)
    deallocate (p, uf, vf)
    deallocate (uc, vc)
    deallocate (pprime, uf_old, vf_old)
    deallocate (ap_p, ae_p, aw_p, an_p, as_p, b_p)
    deallocate (ap_u, ae_u, aw_u, an_u, as_u, b_u, d_u)
    deallocate (ap_v, ae_v, aw_v, an_v, as_v, b_v, d_v)
    
  end subroutine deallocate_memory
  
  subroutine initialize()
    implicit none
    integer :: i, j 
    
    ! generate grid
    dx = Lx/float(nx)
    dy = Ly/float(ny)
    
    do j = 0, ny+1
       do i = 1, nx
          xgrid(i, j) = (i-1+0.5)*dx
       end do
    end do
    
    do j = 1, ny
       do i = 0, nx+1
          ygrid(i, j) = (j-1+0.5)*dy
       end do
    end do
    
    xgrid(0, :) = 0.0d0
    !xgrid(nx+1, :) = (nx-1+0.5)*dx
    xgrid(nx+1, :) = Lx
    
    ygrid(:, 0) = 0.0d0
    !ygrid(:, ny+1) = (ny-1+0.5)*dy
    ygrid(:, ny+1) = Ly
    
    ! initialize variables
    p = 0.0d0
    uf = 0.0d0
    vf = 0.0d0
    
    ! set bc - top boundary
    uf(:, ny+1) = ulid
    vf(:, ny)   = 0.0d0
    
    ! set bc - bottom boundary
    uf(:, 0)    = 0.0d0
    vf(:, 0)    = 0.0d0

    ! set bc - left boundary
    uf(0, :)    = 0.0d0
    vf(0, :)    = 0.0d0

    ! set bc - right boundary
    uf(nx, :)   = 0.0d0
    vf(nx+1, :) = 0.0d0

    uf_old = uf
    vf_old = vf

  end subroutine initialize

  subroutine fill_u_momentum_coeffs()
    implicit none
    integer :: i, j, k
    double precision :: De, Dw, Dn, Ds
    double precision :: Fe, Fw, Fn, Fs

    ! u-face storage
    ! i = 0 is left-boundary
    ! i = nx is right-boundary 
    ! i = 1 to nx-1 are interior u-cell values in x-direction
    ! j = 1 to ny are interior u-cell values in y-direction
    
    do j = 1, ny
       do i = 1, nx-1
          
          De = mu*dy/dx
          Dw = mu*dy/dx
          Dn = mu*dx/dy
          Ds = mu*dx/dy
          
          Fe = rho*0.50d0*(uf_old(i, j)+uf_old(i+1, j))*dy
          Fw = rho*0.50d0*(uf_old(i-1, j)+uf_old(i, j))*dy
          Fn = rho*0.50d0*(vf_old(i, j)+vf_old(i+1, j))*dx
          Fs = rho*0.50d0*(vf_old(i, j-1)+vf_old(i+1, j-1))*dx

          if (j == 1) Ds = 2.0d0*Ds
          if (j == ny) Dn = 2.0d0*Dn

          ae_u(i, j) = De + max(-Fe, 0.0d0)
          aw_u(i, j) = Dw + max( Fw, 0.0d0)
          an_u(i, j) = Dn + max(-Fn, 0.0d0)
          as_u(i, j) = Ds + max( Fs, 0.0d0)
          
          ap_u(i, j) = ae_u(i, j) + aw_u(i, j) + &
               &       an_u(i, j) + as_u(i, j) + &
               &       Fe - Fw + Fn - Fs

          b_u(i, j) = dy*(p(i, j) - p(i+1, j))
          
          d_u(i, j) = alpha_u*dy/ap_u(i, j)
          
       end do
    end do

  end subroutine fill_u_momentum_coeffs

  subroutine fill_v_momentum_coeffs()
    implicit none
    integer :: i, j
    double precision :: De, Dw, Dn, Ds
    double precision :: Fe, Fw, Fn, Fs

    ! v-face storage              
    ! i = 0 is out-side left-boundary (referred to as left-boundary)
    ! i = nx+1 is out-side right-boundary (referred to as right-boundary value)
    ! j = 0 is bottom boundary 
    ! j = ny is top boundary 
    ! i = 1 to nx are interior values 
    ! j = 1 to ny-1 are interior values 

    do j = 1, ny-1
       do i = 1, nx

          De = mu*dy/dx
          Dw = mu*dy/dx
          Dn = mu*dx/dy
          Ds = mu*dx/dy

          Fe = rho*0.50d0*(uf_old(i, j) + uf_old(i, j+1))*dy
          Fw = rho*0.50d0*(uf_old(i-1, j) + uf_old(i-1, j+1))*dy
          Fn = rho*0.50d0*(vf_old(i, j)+vf_old(i, j+1))*dx
          Fs = rho*0.50d0*(vf_old(i, j-1)+vf_old(i, j))*dx

          if (i == 1) Dw = 2.0d0*Dw
          if (i == nx) De = 2.0d0*De
          
          ae_v(i, j) = De + max(-Fe, 0.0d0)
          aw_v(i, j) = Dw + max( Fw, 0.0d0)
          an_v(i, j) = Dn + max(-Fn, 0.0d0)
          as_v(i, j) = Ds + max( Fs, 0.0d0)

          ap_v(i, j) = ae_v(i, j) + aw_v(i, j) + &
               &       an_v(i, j) + as_v(i, j) + &
               &       Fe - Fw + Fn - Fs

          b_v(i, j) = dx*(p(i, j) - p(i, j+1))

          d_v(i, j) = dx*alpha_v/ap_v(i, j)

       end do
    end do

  end subroutine fill_v_momentum_coeffs
  
  subroutine fill_pressure_correction_coeffs()
    implicit none
    integer :: i, j
    double precision :: error 
    double precision :: Fe, Fw, Fn, Fs

    do j = 1, ny
       do i = 1, nx

          Fe = rho*dy*uf(i, j)
          Fw = rho*dy*uf(i-1, j)
          Fn = rho*dx*vf(i, j)
          Fs = rho*dx*vf(i, j-1)

          ae_p(i, j) = rho*dy*d_u(i, j)
          aw_p(i, j) = rho*dy*d_u(i-1, j)
          an_p(i, j) = rho*dx*d_v(i, j)
          as_p(i, j) = rho*dx*d_v(i, j-1)

          b_p(i, j) = -1.0d0*(Fe - Fw + Fn - Fs)

       end do
    end do

    error = sum(b_p**2.0)
    error = sqrt(error)

    write(*, '(a, f15.10)') 'mass imbalance before solving PCE = ', error

    ! update coefficients along boundaries
    ae_p(nx, :) = 0.0d0
    aw_p(1, :)  = 0.0d0
    an_p(:, ny) = 0.0d0
    as_p(:, 1)  = 0.0d0

    ap_p = ae_p + aw_p + an_p + as_p
    ap_p(1, 1) = 1.0d40         ! setting reference cell pressure

  end subroutine fill_pressure_correction_coeffs

  subroutine solve_u_momentum_equation_gs()
    implicit none
    integer :: i, j, iter
    double precision :: error
    
    do iter = 1, u_iter
       do j = 1, ny
          do i = 1, nx-1
             
             uf(i, j) = (ae_u(i, j)*uf(i+1, j) + aw_u(i, j)*uf(i-1, j) + &
                  &     an_u(i, j)*uf(i, j+1) + as_u(i, j)*uf(i, j-1) + b_u(i, j) + &
                  &     ((1.0d0 - alpha_u)/alpha_u)*uf_old(i, j)*ap_u(i, j))*(alpha_u/ap_u(i, j))

          end do
       end do
    end do

    ! compute u-error
    error = 0.0d0
    do j = 1, ny
       do i = 1, nx-1
          
          error = error + (uf(i, j) - uf_old(i, j))**2.0
          
       end do
    end do

    error = sqrt(error)
    
    write(*, '(a, f15.10)') 'U-momentum error = ', error

  end subroutine solve_u_momentum_equation_gs

  subroutine solve_v_momentum_equation_gs()
    implicit none
    integer :: i, j, iter
    double precision :: error

    do iter = 1, v_iter       
       do j = 1, ny-1
          do i = 1, nx 
             
             vf(i, j) = (ae_v(i, j)*vf(i+1, j) + aw_v(i, j)*vf(i-1, j) + &
                  &     an_v(i, j)*vf(i, j+1) + as_v(i, j)*vf(i, j-1) + b_v(i, j) + &
                  &     ((1.0d0 - alpha_v)/alpha_v)*vf_old(i, j)*ap_v(i, j))*(alpha_v/ap_v(i, j))

          end do
       end do
    end do
    
    error = 0.0d0
    do j = 1, ny-1
       do i = 1, nx
          
          error = error + (vf(i, j) - vf_old(i, j))**2.0
          
       end do
    end do

    error = sqrt(error)
    
    write(*, '(a, f15.10)') 'V-momentum error = ', error

  end subroutine solve_v_momentum_equation_gs

  subroutine solve_pressure_correction_equation_gs()
    implicit none
    integer :: i, j, iter
    double precision :: error

    ! important to initialize pprime every time to take care of wall BCs
    pprime = 0.0d0

    do iter = 1, p_iter
       do j = 1, ny
          do i = 1, nx 
             
             pprime(i, j) = (ae_p(i, j)*pprime(i+1, j) + aw_p(i, j)*pprime(i-1, j) + &
                  &     an_p(i, j)*pprime(i, j+1) + as_p(i, j)*pprime(i, j-1) + b_p(i, j)) &
                  &     /ap_p(i, j)
          end do
       end do
    end do

  end subroutine solve_pressure_correction_equation_gs

  subroutine correct_pressure_and_velocity()
    implicit none
    integer :: i, j
    double precision :: bstar

    ! correct u-velocity
    do j = 1, ny
       do i = 1, nx-1
          uf(i, j) = uf(i, j) + d_u(i, j)*(pprime(i, j) - pprime(i+1, j))
       end do
    end do

    ! correct v-velocity
    do j = 1, ny-1
       do i = 1, nx
          vf(i, j) = vf(i, j) + d_v(i, j)*(pprime(i, j) - pprime(i, j+1))
       end do
    end do

    ! correct pressure
    do j = 1, ny
       do i = 1, nx
          p(i, j) = p(i, j) + (pprime(i, j))*alpha_p
       end do
    end do

    bstar = 0.0d0
    do j = 1, ny
       do i = 1, nx
          bstar = bstar + rho*uf(i, j)*dy- rho*uf(i-1, j)*dy + rho*vf(i, j)*dx - rho*vf(i, j-1)*dx
       end do
    end do

    uf_old = uf
    vf_old = vf

    write(*, '(a, f15.10)') 'mass imbalance after solving PCE = ', bstar

  end subroutine correct_pressure_and_velocity

  subroutine print_diagnostics()
    implicit none
    !    integer :: i, j

    write(*, '(3(f15.10, 5x))') uf(nx/2, ny/2), vf(nx/2, ny/2), p(nx/2, ny/2)

  end subroutine print_diagnostics

  subroutine calculate_residuals()
    implicit none
    integer :: i, j 
    double precision :: Fe, Fw, Fn, Fs

    residual = 0.0d0

    ! calculate u-residual
    u_residual = 0.0d0

    do j = 1, ny
       do i = 1, nx-1 

          u_residual = u_residual -  ap_u(i, j)*uf(i, j)/alpha_u + &
               (ae_u(i, j)*uf(i+1, j) + aw_u(i, j)*uf(i-1, j) + &
               &     an_u(i, j)*uf(i, j+1) + as_u(i, j)*uf(i, j-1) + b_u(i, j) + &
               &     ((1.0d0 - alpha_u)/alpha_u)*uf(i, j)*ap_u(i, j))

       end do
    end do

    v_residual = 0.0d0

    do j = 1, ny-1
       do i = 1, nx 

          v_residual = v_residual - ap_v(i, j)*vf(i, j)/alpha_v + &
               &     (ae_v(i, j)*vf(i+1, j) + aw_v(i, j)*vf(i-1, j) + &
               &     an_v(i, j)*vf(i, j+1) + as_v(i, j)*vf(i, j-1) + b_v(i, j) + &
               &     ((1.0d0 - alpha_v)/alpha_v)*vf(i, j)*ap_v(i, j))

       end do
    end do

    p_residual = 0.0d0

    do j = 1, ny
       do i = 1, nx

          Fe = rho*dy*uf(i, j)
          Fw = rho*dy*uf(i-1, j)
          Fn = rho*dx*vf(i, j)
          Fs = rho*dx*vf(i, j-1)

          p_residual = p_residual + Fe - Fw + Fn - Fs

       end do
    end do

    residual = abs(u_residual) + abs(v_residual) + abs(p_residual)

  end subroutine calculate_residuals

  subroutine write_output_file()
    implicit none
    integer :: i, j

    ! interpolate velocity to cell centers
    do j = 1, ny
       do i = 1, nx

          uc(i, j) = 0.50d0*(uf(i-1, j) + uf(i, j))
          vc(i, j) = 0.50d0*(vf(i, j) + vf(i, j-1))
          
       end do
    end do

    uc(0, 0:ny+1) = uf(0, 0:ny+1)
    uc(nx+1, 0:ny+1) = uf(nx, 0:ny+1)
    uc(0:nx, 0) = uf(0:nx, 0)
    uc(0:nx, ny+1) = uf(0:nx, ny+1)
    
    vc(0:nx+1, 0) = vf(0:nx+1, 0)
    vc(0:nx+1, ny+1) = vf(0:nx+1, ny)
    vc(0, 0:ny) = vf(0, 0:ny)
    vc(nx+1, 0:ny) = vf(nx+1, 0:ny)

    ! update boundary pressure with homogeneous Neumann BC
    do i = 0, nx+1
       p(i, 0) = p(i, 1)
       p(i, ny+1) = p(i, ny)
    end do

    do j = 0, ny+1
       p(0, j) = p(1, j)
       p(nx+1, j) = p(nx, j)
    end do

    open(unit = 20, file = 'plot.dat', status = 'replace')
    WRITE(20,*) 'VARIABLES =X, Y, U, V, P'
    WRITE(20,*) 'ZONE ','I=',nx+2,'J=',ny+2,',','F=BLOCK'

    ! output x horizontally
    DO j = 0, ny+1
       WRITE(20,*) (xgrid(i, j),i=0,nx+1)
    END DO
    ! output y vertically
    DO j = 0, ny+1
       WRITE(20,*)(ygrid(i, j),i = 0,nx+1)
    END DO
    ! output u for each node
    DO j = 0, ny+1
       WRITE(20,*)(uc(i,j),i=0,nx+1)
    END DO
    ! output v for each node
    DO j= 0, ny+1
       WRITE(20,*)(vc(i,j),i=0,nx+1)
    END DO
    ! output p for each node
    DO j= 0, ny+1
       WRITE(20,*)(p(i,j),i=0,nx+1)
    END DO
    close(20)

    write(*, *) 'Wrote output file'

  end subroutine write_output_file

end module data

program main 
  use data
  implicit none
  integer :: iter, i, j
  double precision :: tol = 1.0d-6

  residual = 1.0d0

  call allocate_memory()

  call initialize()

  iter = 1

  do while (residual > tol)

     write(*, '(a, i7)') "SIMPLE iteration = ", iter
     
     call fill_u_momentum_coeffs()
     call fill_v_momentum_coeffs()

     call solve_u_momentum_equation_gs()
     call solve_v_momentum_equation_gs()

     call fill_pressure_correction_coeffs()
     call solve_pressure_correction_equation_gs()

     call correct_pressure_and_velocity()

     call calculate_residuals()

     write(*, '(a, f15.10)') 'Global residual = ', residual

     write(*, '(a)') 'Printing diagnostics'

     call print_diagnostics()
     
     write(70, *) iter, residual

     ! if (mod(iter,50) == 0) then
     !    call write_output_file()
     ! end if

     iter = iter + 1

     write(*, *) '-----------------------------------------------------------------'
     
  end do

  write(*, '(i7, 4(f15.10, 5x))') iter, uf(nx-2, ny-2), vf(nx-2, ny-2), p(nx-2, ny-2), residual

  call write_output_file()
  
  do j = 1, ny+1
     write(50, *) ygrid(nx/2, j), uf(nx/2, j)
  end do
  
  do i = 1, nx
     write(60, *) xgrid(i, ny/2), vf(i, ny/2)
  end do

  ! do j = 1, ny
  !    do i = 1, nx
  !       write(70, *) xgrid(i, j), ygrid(i, j), uf(i, j), vf(i, j), p(i, j)
  !    end do
  ! end do

  call deallocate_memory()

end program main
