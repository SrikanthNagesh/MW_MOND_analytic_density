!!!!!!!!!!! ----- This file is added by STN ---------!!!!!!

module mond_disk_parameters
  use amr_parameters

  implicit none

  ! MOND constants
  real(dp) :: G_mond  = 1.0d0
  real(dp) :: a0_mond = 3.703d0

  ! Galaxy parameters
  real(dp) :: mass_disk = 9.15d10

  real(dp) :: f_inner = 0.8236d0
  real(dp) :: f_outer = 0.1764d0

  real(dp) :: rd_inner = 1.29d0
  real(dp) :: rd_outer = 4.20d0

  real(dp) :: rmax

end module mond_disk_parameters
