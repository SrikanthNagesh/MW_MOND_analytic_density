!########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine rho_ana(x,d,dx,ncell)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  use mond_parameters
  use mond_disk_density
  implicit none
  integer ::ncell                         ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector)       ::d ! Density
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates analytical Poisson source term.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! d(i) is the density field in user units.
  !================================================================
  integer::i
  real(dp)::dmass,emass,xmass,ymass,zmass,rr,rx,ry,rz,dd,rcut
  real(dp):: Mpl, Pi

  Pi = acos(-1.0d0)
  
  emass=2.*boxlen*0.5d0**nlevelmax
  xmass=boxlen/2.0
  ymass=boxlen/2.0
  zmass=boxlen/2.0
  !rcut =boxlen/2.0*0.75

  !rmax=gravity_params(1)
  !dmass=1.0/emass/(1.0+emass)**2

  do i=1,ncell
     rx=x(i,1)-xmass
     ry=x(i,2)-ymass
     rz=x(i,3)-zmass
     rr=sqrt(rx**2+ry**2+rz**2)
     if (rr .le. rmax) then
        call interpolate_density(rr,d(i))
     !if (rr .le. 10.0) then
     !   write(*,*)"The density at given r is", rr, d(i)
     !end if 
     else 
        d(i) = 1.0d-5
     end if 
  end do


end subroutine rho_ana
