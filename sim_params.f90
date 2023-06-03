module sim_params
  ! defines variables for the simulation
  use decomp_2d, only : mytype

  implicit none

  integer, parameter :: nx=101, ny=101, nz=101
  real(mytype), parameter :: xlx=10.0, yly=10.0, zlz=10.0
  real(mytype) :: dx = (xlx - 0)/(nx - 1)
  real(mytype) :: dy = (yly - 0)/(ny - 1)
  real(mytype) :: dz = (zlz - 0)/(nz - 1)
  real(mytype), parameter :: dt=0.00001
  integer :: itime
  integer, parameter :: prow=2, pcol=4
end module sim_params
