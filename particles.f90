module particles
! 
! Author : Nikhil Jayswal, Indian Institute of Science
!
! This module implements particle tracking using the 2DECOMP&FFT library
! It assumes an x-pencil orientation with default origin + coordinate system
!
  use MPI
  use decomp_2d 

  ! only for tracking_v2
  use sim_params

  implicit none
  
  ! variables required by all subroutines
  integer, parameter :: nps = 10000                     ! number of particles
  real(mytype), dimension(nps, 3), private :: cmatrix  ! particle coordinate matrix
  logical, dimension(nps), private :: ppmap = .false.  ! boolean particle-process map
  integer :: local_pcount                              ! local particle count
  integer, dimension(prow*pcol, 2) :: domain_table     ! global domain index table

contains

  subroutine init_particle_tracking
  ! pre-tracking setup

    use MPI
    implicit none
    
    integer, dimension(2) :: domain_limits
    integer :: i, counter, ierror
    integer, dimension(prow*pcol*2) :: domain_array
    
    ! create domain table
    ! row i = (xstart(2), xstart(3)) of process with rank = i
    domain_limits(1) = xstart(2)
    domain_limits(2) = xstart(3)

    call MPI_ALLGATHER(domain_limits, 2, MPI_INT, domain_array, 2, &
              MPI_INT, MPI_COMM_WORLD, ierror)

    counter = 1
    do i = 1, nproc
      domain_table(i, :) = (/domain_array(counter), domain_array(counter + 1)/)
      counter = counter + 2
    end do

  end subroutine init_particle_tracking

  subroutine load_cmatrix
  ! loads presaved coordinate matrix
    implicit none
  end subroutine load_cmatrix
  
  subroutine inject_particles
  ! injects particles into the flow field 
  ! flag particles which are in local domain 
    implicit none

    integer :: i, n
    !integer :: ccount, count_rate, count_max
    integer, allocatable, dimension(:) :: seed
    real(mytype) :: y0, y1, z0, z1

    !call system_clock(ccount, count_rate, count_max) 
    call random_seed(size=n)
    allocate(seed(n))
    seed = 16961
    call random_seed(put=seed)

    do i = 1, nps
      !call random_seed(size=n)
      !allocate(seed(n))
      !seed = 16961
      !call random_seed(put=seed)
      call random_number(cmatrix(i, 1))
      cmatrix(i, 1) = xlx * cmatrix(i, 1)
      cmatrix(i, 2) = yly - 0.000001      ! particles in grid interior
     
      !call random_seed(put=seed)
      call random_number(cmatrix(i, 3))
      cmatrix(i, 3) = zlz*cmatrix(i, 3)
      !deallocate(seed)
      
      !if (nrank == 3) then
      !  print *, cmatrix(i, 1), cmatrix(i, 2), cmatrix(i, 3)
      !endif
 
      ! create/update local particle map after injection
      y0 = (xstart(2) - 1)*dy
      y1 = xend(2)*dy
      z0 = (xstart(3) - 1)*dz
      z1 = xend(3)*dz
      
      if ((cmatrix(i, 2) .ge. y0) .and. (cmatrix(i, 2) .le. y1)) then
          if((cmatrix(i, 3) .ge. z0) .and. (cmatrix(i, 3) .le. z1)) then
              ppmap(i) = .true.
          end if
      else
        ppmap(i) = .false.
      end if

    end do

      ! set local particle count
      local_pcount = count(ppmap)

    deallocate(seed)
  end subroutine inject_particles

  subroutine write_particles(iter)
    use MPI

    implicit none

    integer :: ierror, fh, i, j, iter, type_size
    character(len=80) :: fname
    real(mytype), allocatable, dimension(:) :: write_cmatrix
    integer, dimension(nproc) :: pcount_array
    integer(kind=MPI_OFFSET_KIND) :: fdisp

    ! construct file name
    write(fname, '(a, i0)') 'cmatrix', iter

    if (local_pcount .gt. 0) then
      allocate(write_cmatrix(local_pcount*4))    ! allocate memory to write data
    else
      allocate(write_cmatrix(0))
    end if
    
    ! create particle location data array
    j = 1
    do i = 1, nps
      if (ppmap(i)) then
        write_cmatrix(j) = i*1.0
        write_cmatrix(j+1) = cmatrix(i, 1)
        write_cmatrix(j+2) = cmatrix(i, 2)
        write_cmatrix(j+3) = cmatrix(i, 3)
        j = j + 4
      end if
    end do

      !do j = 1, local_pcount
      !  print *, write_cmatrix(4*j - 3 : 4*j)
      !end do

    ! open file for writing
    call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(fname), MPI_MODE_CREATE+MPI_MODE_WRONLY, &
                                MPI_INFO_NULL, fh, ierror)
  
    ! get particle counts from each proc
    call MPI_ALLGATHER(local_pcount, 1, MPI_INT, pcount_array, 1, MPI_INT, &
                                MPI_COMM_WORLD, ierror)

!    if (nrank == 0) then
!      print *, pcount_array
!    end if
    ! compute displacements
    fdisp = 0
    if (nrank .gt. 0) then
      do i = 1, nrank
          fdisp = fdisp + pcount_array(i)
      end do
    end if
    call MPI_TYPE_SIZE(MPI_DOUBLE, type_size, ierror)
    fdisp = fdisp*4*type_size

    !print *, 'nrank = ', nrank, 'pcount = ', local_pcount, 'offset = ', fdisp
    ! write particle location data to file
    !call MPI_FILE_SEEK(fh, fdisp*4, MPI_SEEK_SET, ierror)
    call MPI_FILE_WRITE_AT(fh, fdisp, write_cmatrix, local_pcount*4, &
                            MPI_DOUBLE, MPI_STATUS_IGNORE, ierror)
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
    
    ! clean up
    deallocate(write_cmatrix)
    call MPI_FILE_CLOSE(fh, ierror)
  end subroutine write_particles 

  subroutine trilinear_interp(u_x, u_y, u_z, p_coords, vel_x, vel_y, vel_z)
    use decomp_2d

    implicit none

    real(mytype), dimension(1:xsize(1), 0:xsize(2)+1, 0:xsize(3)+1) :: u_x,u_y,u_z ! halo fields
    real(mytype), dimension(3), intent(in) :: p_coords
    integer :: i1, j1, k1
    real(mytype) :: x0, x1, y0, y1, z0, z1
    real(mytype) :: xd, yd, zd
    real(mytype) :: c00, c01, c10, c11
    real(mytype) :: c0, c1, c_x, c_y, c_z
    real(mytype), intent(out) :: vel_x, vel_y, vel_z
   
    c_x = p_coords(1)
    c_y = p_coords(2)
    c_z = p_coords(3)
    
    i1 = floor(c_x/dx)
    j1 = floor(c_y/dy)
    k1 = floor(c_z/dz)

    x1 = (i1 + 1)*dx
    x0 = i1*dx
    y1 = (j1 + 1)*dy
    y0 = j1*dy
    z1 = (k1 + 1)*dz
    z0 = k1*dz
    
    xd = (c_x - x0)/(x1 - x0)
    yd = (c_y - y0)/(y1 - y0)
    zd = (c_z - z0)/(z1 - z0)
    
    ! local index - left one - min. value = 1
    i1 = floor(c_x/dx) - (xstart(1) - 1) + 1
    j1 = floor(c_y/dy) - (xstart(2) - 1) + 1
    k1 = floor(c_z/dz) - (xstart(3) - 1) + 1

    c00 = u_x(i1, j1, k1)*(1.0 - xd) + u_x(i1 + 1, j1, k1)*xd
    c01 = u_x(i1, j1, k1 + 1)*(1.0 - xd) + u_x(i1 + 1, j1, k1 + 1)*xd
    c10 = u_x(i1, j1 + 1, k1)*(1.0 - xd) + u_x(i1 + 1, j1 + 1, k1)*xd
    c11 = u_x(i1, j1 + 1, k1 + 1)*(1.0 - xd) + u_x(i1 + 1, j1 + 1, k1 + 1)*xd

    c0 = c00*(1.0 - yd) + c10*yd
    c1 = c01*(1.0 - yd) + c11*yd
   
    vel_x = c0*(1.0 - zd) + c1*zd

    c00 = u_y(i1, j1, k1)*(1.0 - xd) + u_y(i1 + 1, j1, k1)*xd
    c01 = u_y(i1, j1, k1 + 1)*(1.0 - xd) + u_y(i1 + 1, j1, k1 + 1)*xd
    c10 = u_y(i1, j1 + 1, k1)*(1.0 - xd) + u_y(i1 + 1, j1 + 1, k1)*xd
    c11 = u_y(i1, j1 + 1, k1 + 1)*(1.0 - xd) + u_y(i1 + 1, j1 + 1, k1 + 1)*xd

    c0 = c00*(1.0 - yd) + c10*yd
    c1 = c01*(1.0 - yd) + c11*yd
   
    vel_y = c0*(1.0 - zd) + c1*zd

    c00 = u_z(i1, j1, k1)*(1.0 - xd) + u_z(i1 + 1, j1, k1)*xd
    c01 = u_z(i1, j1, k1 + 1)*(1.0 - xd) + u_z(i1 + 1, j1, k1 + 1)*xd
    c10 = u_z(i1, j1 + 1, k1)*(1.0 - xd) + u_z(i1 + 1, j1 + 1, k1)*xd
    c11 = u_z(i1, j1 + 1, k1 + 1)*(1.0 - xd) + u_z(i1 + 1, j1 + 1, k1 + 1)*xd

    c0 = c00*(1.0 - yd) + c10*yd
    c1 = c01*(1.0 - yd) + c11*yd
   
    vel_z = c0*(1.0 - zd) + c1*zd
  end subroutine trilinear_interp

  subroutine track_particles(u_x, u_y, u_z)
  ! marches particles forward one time-step using RK2
  ! transfers particles crossing domains

    implicit none

    integer :: i, j, k, j1, k1, counter
    integer, dimension(2) :: row
    real(mytype) :: p_coords(3)
    real(mytype) :: vel_x, vel_y, vel_z
    real(mytype) :: c_x, c_y, c_z
    real(mytype) :: y0, y1, z0, z1
    real(mytype), dimension(1:xsize(1),0:xsize(2)+1,0:xsize(3)+1) :: u_x, u_y, u_z ! halo fields
    integer :: transfer_count_sum = 0
    integer, dimension(nproc) :: transfer_count_array
    integer :: ierror
    integer :: transfer_count = 0                          
    integer, dimension(nproc) :: displs
    real(mytype), allocatable, dimension(:) :: all_transfer_data
    integer :: c, p_id
    real(mytype), allocatable, dimension(:) :: transfer_data_array
 
    allocate(transfer_data_array(5*count(ppmap)))
    do i = 1, nps
      if (ppmap(i)) then
        ! march particle forward in time
        p_coords =  cmatrix(i, :)
          
        ! RK2 Step 1
        call trilinear_interp(u_x, u_y, u_z, p_coords, vel_x, vel_y, vel_z)
        p_coords = p_coords + (dt/2)*(/vel_x, vel_y, vel_z/)
      
        ! check if out of domain
        c_x = p_coords(1)
        c_y = p_coords(2)
        c_z = p_coords(3)

        if ((c_x > xlx) .or. (c_y > yly) .or. (c_z > zlz)) then
          ppmap(i) = .false.
          cycle
        end if
       
        if ((c_x < 0) .or. (c_y < 0) .or. (c_z < 0)) then
          ppmap(i) = .false.
          cycle
        end if


        ! RK2 Step 2
        call trilinear_interp(u_x, u_y, u_z, p_coords, vel_x, vel_y, vel_z)
        p_coords = p_coords + (dt/2)*(/vel_x, vel_y, vel_z/)
       
        ! update cmatrix
        cmatrix(i, :) = p_coords

        ! flag for transfer if required
        c_x = p_coords(1)
        c_y = p_coords(2)
        c_z = p_coords(3)

        ! rule out particles that have wandered off global domain limits
        if ((c_x > xlx) .or. (c_y > yly) .or. (c_z > zlz)) then
          ppmap(i) = .false.
          cycle
        end if
       
        if ((c_x < 0) .or. (c_y < 0) .or. (c_z < 0)) then
          ppmap(i) = .false.
          cycle
        end if

        y1 = dy*xend(2) ! local domain extended by one halo cell in one dir'n
        z1 = dz*xend(3) ! no boundary issues, no particles out of global domain
        y0 = dy*(xstart(2) - 1)
        z0 = dz*(xstart(3) - 1)

        if ((c_y .gt. y1) .or. (c_z .gt. z1) .or. (c_y .lt. y0) .or. (c_z .lt. z0)) then
          ppmap(i) = .false.
          ! find new proc - whose xstart(2, 3) = nearest(floor) y,z node
          j1 = floor(c_y/dy) + 1 
          k1 = floor(c_z/dz) + 1
          j = maxval(domain_table(:, 1), mask = domain_table(:, 1) .le. j1)
          k = maxval(domain_table(:, 2), mask = domain_table(:, 2) .le. k1)
          row = (/j, k/)
          do j = 1, nproc
            if (all(domain_table(j, :) == row)) then
              k = j
              exit
            end if
          end do

          ! write particle transfer data
          transfer_count = transfer_count + 1
          counter = (transfer_count - 1)*5 + 1
          transfer_data_array(counter) = (k - 1)*1.0
          transfer_data_array(counter + 1) = i*1.0
          transfer_data_array(counter + 2) = c_x
          transfer_data_array(counter + 3) = c_y
          transfer_data_array(counter + 4) = c_z 
        end if
      end if
    end do

    ! compute total number of particles to be transferred
    call MPI_ALLREDUCE(transfer_count, transfer_count_sum, 1, MPI_INT, MPI_SUM, &
                                            MPI_COMM_WORLD, ierror)
    ! can this be changed to a more faster check? min of all > 0?
 
    if (transfer_count_sum > 0) then ! is this necessary? how much overhead when 
                                     ! transfer count = 0?
      ! get individual counts of particles to be transferred
      call MPI_ALLGATHER(transfer_count, 1, MPI_INT, transfer_count_array, 1, &
                                      MPI_INT, MPI_COMM_WORLD, ierror)

      allocate(all_transfer_data(5*transfer_count_sum))
      all_transfer_data = 0.0
      ! compute displacements
      displs(1) = 0
      do i = 2, nproc
        displs(i) = displs(i-1) + transfer_count_array(i-1)
      end do
  
      ! gather all transfer data
      call MPI_ALLGATHERV(transfer_data_array, 5*transfer_count, MPI_DOUBLE, &
            all_transfer_data, 5*transfer_count_array, 5*displs, MPI_DOUBLE, &
                  MPI_COMM_WORLD, ierror)

      ! check for newly assigned particles and copy data
      do i = 1, transfer_count_sum
        counter = (i-1)*5 + 1
        c = int(all_transfer_data(counter))
        if (c == nrank) then
          p_id = int(all_transfer_data(counter + 1))
          cmatrix(p_id, 1) = all_transfer_data(counter + 2)
          cmatrix(p_id, 2) = all_transfer_data(counter + 3)
          cmatrix(p_id, 3) = all_transfer_data(counter + 4)
          ppmap(p_id) = .true.
        end if
      end do

      ! reset transfer count
      transfer_count  = 0

      ! clean up
      deallocate(all_transfer_data)
      deallocate(transfer_data_array)
    end if

    ! update particle count
    local_pcount = count(ppmap)
  end subroutine track_particles

end module
